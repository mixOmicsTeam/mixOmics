####################################################################################
## ---------- internal
.spca <- function(X,
                  ncomp = 2,
                  center = TRUE,
                  keepX = NULL,
                  scale = TRUE,
                  max.iter = 500,
                  tol = 1e-06,
                  logratio = c('none', 'CLR'),
                  multilevel = NULL) {

  ## ----------------------------------- checks
  logratio <- .matchArg(logratio)
  mcd <- mget(names(formals()),sys.frame(sys.nframe())) ## match.call and defaults
  err = tryCatch(mcd, error = function(e) e) ## see if arguments can be evaluated
  if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)

  ## check pca entries for match.call and defaults and do necessary adjustments to
  ## arguments in this environement from within function
  .pcaEntryChecker(mcd, check.keepX = TRUE, check.center = TRUE, check.NA = FALSE)

  ## ----------------------------------- logratio tran.
  X = logratio.transfo(X = X, logratio = logratio, offset = 0)#if(logratio == "ILR") {ilr.offset} else {0})

  ## ----------------------------------- multilevel
  if (!is.null(multilevel))
  {
    # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
    multilevel = data.frame(multilevel)

    if ((nrow(X) != nrow(multilevel)))
      stop("unequal number of rows in 'X' and 'multilevel'.")

    if (ncol(multilevel) != 1)
      stop("'multilevel' should have a single column for the repeated measurements.")

    multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements

    Xw = withinVariation(X, design = multilevel)
    X = Xw
  }

  #--scaling the data--#
  X=scale(X,center=center,scale=scale)
  cen = attr(X, "scaled:center")
  sc = attr(X, "scaled:scale")
  if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance.")

  #--initialization--#
  X=as.matrix(X)
  X.temp=as.matrix(X)
  n=nrow(X)
  p=ncol(X)

  # put a names on the rows and columns
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
  colnames(X) = X.names

  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) ind.names = 1:nrow(X)
  rownames(X) = ind.names

  vect.varX=vector(length=ncomp)
  names(vect.varX) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)

  vect.iter=vector(length=ncomp)
  names(vect.iter) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)

  vect.keepX = keepX
  names(vect.keepX) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)

  # KA: to add if biplot function (but to be fixed!)
  #sdev = vector(length = ncomp)

  mat.u=matrix(nrow=n, ncol=ncomp)
  mat.v=matrix(nrow=p, ncol=ncomp)
  colnames(mat.u)=paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
  colnames(mat.v)=paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
  rownames(mat.v)=colnames(X)

  #--loop on h--#
  for(h in 1:ncomp){

    #--computing the SVD--#
    svd.X = svd(X.temp, nu = 1, nv = 1)
    u = svd.X$u[,1]
    loadings = svd.X$v[,1]#svd.X$d[1]*svd.X$v[,1]
    v.stab = FALSE
    u.stab = FALSE
    iter = 0

    #--computing nx(degree of sparsity)--#
    nx = p-keepX[h]
    #vect.keepX[h] = keepX[h]

    u.old = u
    loadings.old = loadings
    #--iterations on v and u--#
    repeat{

      iter = iter +1

      loadings = t(X.temp) %*% u

      #--penalisation on loading vectors--#
      if (nx != 0) {
        absa = abs(loadings)
        if(any(rank(absa, ties.method = "max") <= nx)) {
          loadings = ifelse(rank(absa, ties.method = "max") <= nx, 0, sign(loadings) * (absa - max(absa[rank(absa, ties.method = "max") <= nx])))
        }
      }
      loadings = loadings / drop(sqrt(crossprod(loadings)))

      u = as.vector(X.temp %*% loadings)
      #u = u/sqrt(drop(crossprod(u))) # no normalisation on purpose: to get the same $x as in pca when no keepX.

      #--checking convergence--#
      if(crossprod(u-u.old)<tol){break}
      if(crossprod(loadings-loadings.old)<tol){break}

      if (iter >= max.iter)
      {
        warning(paste("Maximum number of iterations reached for the component", h),call. = FALSE)
        break
      }

      u.old = u
      loadings.old = loadings

    }
    #v.final = v.new/sqrt(drop(crossprod(v.new)))

    #--deflation of data--#
    c = crossprod(X.temp, u) / drop(crossprod(u))
    X.temp= X.temp - u %*% t(c)  #svd.X$d[1] * svd.X$u[,1] %*% t(svd.X$v[,1])

    vect.iter[h] = iter
    mat.v[,h] = loadings
    mat.u[,h] = u

    #--calculating adjusted variances explained--#
    X.var = X %*% mat.v[,1:h]%*%solve(t(mat.v[,1:h])%*%mat.v[,1:h])%*%t(mat.v[,1:h])
    vect.varX[h] = sum(X.var^2)

    # KA: to add if biplot function (but to be fixed!)
    #sdev[h] = sqrt(svd.X$d[1])

  }#fin h

  rownames(mat.u) = ind.names
  {
    ## keeping this for the exported generic's benefit, but is not necessary as far as methods are concerned
    mcr = match.call() ## match call for returning
    mcr[[1]] = as.name('spca')
    mcr[-1L] <- lapply(mcr[-1L], eval.parent)
  }

  result = list(
    call = mcr,
    X = X,
    ncomp = ncomp,
    #sdev = sdev,  # KA: to add if biplot function (but to be fixed!)
    #center = center, # KA: to add if biplot function (but to be fixed!)
    #scale = scale,   # KA: to add if biplot function (but to be fixed!)
    varX = vect.varX / sum(X ^ 2),
    keepX = vect.keepX,
    iter = vect.iter,
    rotation = mat.v,
    x = mat.u,
    names = list(X = X.names, sample = ind.names),
    loadings = list(X = mat.v),
    variates = list(X = mat.u)
  )

  class(result) = c("spca", "prcomp", "pca")

  #calcul explained variance
  explX=explained_variance(X,result$variates$X,ncomp)
  result$explained_variance=explX


  return(invisible(result))
}


#' @title Sparse Principal Components Analysis
#'
#' @description Performs a sparse principal components analysis to perform variable
#' selection by using singular value decomposition.
#'
#' The calculation employs singular value decomposition of the (centered and
#' scaled) data matrix and LASSO to generate sparsity on the loading vectors.
#'
#' \code{scale= TRUE} is highly recommended as it will help obtaining
#' orthogonal sparse loading vectors.
#'
#' \code{keepX} is the number of variables to keep in loading vectors. The
#' difference between number of columns of \code{X} and \code{keepX} is the
#' degree of sparsity, which refers to the number of zeros in each loading
#' vector.
#'
#' Note that \code{spca} does not apply to the data matrix with missing values.
#' The biplot function for \code{spca} is not available.
#'
#' According to Filzmoser et al., a ILR log ratio transformation is more
#' appropriate for PCA with compositional data. Both CLR and ILR are valid.
#'
#' Logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#'
#' Logratio can only be applied if the data do not contain any 0 value (for
#' count data, we thus advise the normalise raw data with a 1 offset). For ILR
#' transformation and additional offset might be needed.

## ----------------------------------- Parameters
#' @inheritParams pca
#' @param keepX numeric vector of length ncomp, the number of variables to keep.
#' @param logratio one of ('none','CLR'). Specifies the log ratio
#' transformation to deal with compositional values that may arise from
#' specific normalisation in sequencing data. Default to 'none'.
#' in loading vectors. By default all variables are kept in the model. See
#' details.

## ----------------------------------- Value
#' @return \code{spca} returns a list with class \code{"spca"} containing the
#' following components: \item{ncomp}{the number of components to keep in the
#' calculation.} \item{varX}{the adjusted cumulative percentage of variances
#' explained.} \item{keepX}{the number of variables kept in each loading
#' vector.} \item{iter}{the number of iterations needed to reach convergence
#' for each component.} \item{rotation}{the matrix containing the sparse
#' loading vectors.} \item{x}{the matrix containing the principal components.}
#' Refer to \code{pca} for different type of inputs supported.

## ----------------------------------- Misc
#' @author Ignacio Gonzalez, Kim-Anh Lê Cao, Fangzhou Yao, Al J Abadi
#' @seealso \code{\link{pca}}, \code{\link{ipca}}, \code{\link{selectVar}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}} and http://www.mixOmics.org
#' for more details.
#' @keywords algebra

## ----------------------------------- Examples
#' @example examples/spca-example.R
####################################################################################
## ---------- Generic
#' @param ... Aguments passed to the generic.
#' @export
spca <- function(X=NULL, data=NULL, ncomp=2, keepX=NULL, ...) UseMethod('spca')

####################################################################################
## ---------- Methods

## ----------------------------------- X=matrix
#' @rdname spca
#' @export
spca.default <- .spca

## ----------------------------------- X= assay name from data
#' @importFrom SummarizedExperiment assay
#' @rdname spca
#' @export
spca.character <- function(X=NULL, data=NULL, ncomp=2, keepX=NULL, ...){
  .pcaMethodsHelper(match.call(), fun = 'spca')
}
