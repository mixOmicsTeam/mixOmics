## ----------- Description ----------- 
#' Sparse Independent Principal Component Analysis
#'
#' Performs sparse independent principal component analysis on the given data
#' matrix to enable variable selection.
#'
#' See Details of ipca.
#'
#' Soft thresholding is implemented on the independent loading vectors to
#' obtain sparse loading vectors and enable variable selection.
#'
## ----------- Parameters ----------- 
## keepX
#' @inheritParams spca
## other args
#' @inheritParams ipca
#' 
## ----------- Value ----------- 
#' @return \code{pca} returns a list with class \code{"ipca"} containing the
#' following components: \item{ncomp}{the number of principal components used.}
#' \item{unmixing}{the unmixing matrix of size (ncomp x ncomp)}
#' \item{mixing}{the mixing matrix of size (ncomp x ncomp} \item{X}{the
#' centered data matrix} \item{x}{the principal components (with sparse
#' independent loadings)} \item{loadings}{the sparse independent loading
#' vectors} \item{kurtosis}{the kurtosis measure of the independent loading
#' vectors}
#' 
## ----------- Ref ----------- 
#' @author Fangzhou Yao and Jeff Coquery.
#' @seealso \code{\link{ipca}}, \code{\link{pca}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}} and http://www.mixOmics.org for more details.
#' @references Yao, F., Coquery, J. and Lê Cao, K.-A. (2011) Principal
#' component analysis with independent loadings: a combination of PCA and ICA.
#' (in preparation)
#'
#' A. Hyvarinen and E. Oja (2000) Independent Component Analysis: Algorithms
#' and Applications, \emph{Neural Networks}, \bold{13(4-5)}:411-430
#'
#' J L Marchini, C Heaton and B D Ripley (2010). fastICA: FastICA Algorithms to
#' perform ICA and Projection Pursuit. R package version 1.1-13.
#' @keywords algebra
#' 
## ----------- Examples ----------- 
#' @example examples/sipca-example.R
## setting the document name here so internal would not force the wrong name
#' @name sipca
NULL
## ----------- Internal ----------- 
.sipca <- function(X,
                   ncomp  = 3,
                   keepX = NULL,
                   mode = c("deflation", "parallel"),
                   fun = c("logcosh", "exp"),
                   scale = FALSE,
                   max.iter = 200,
                   tol = 1e-04,
                   w.init = NULL,
                   ret.call=FALSE) {
  mode <- .matchArg(mode)
  fun <- .matchArg(fun)
  
  mcd <- mget(names(formals()),sys.frame(sys.nframe())) ## match.call and defaults
  err = tryCatch(mcd, error = function(e) e) ## see if arguments can be evaluated
  if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)
  
  ## check entries for args and do necessary adjustments to them
  mcd <- .pcaEntryChecker(mcd,fun = 'sipca')
  
  X <- mcd$X
  keepX <- mcd$keepX
  ncomp <- mcd$ncomp
  max.iter <- as.integer(mcd$max.iter)
  rm(mcd)
  
  dim_x <- dim(X)
  d <- dim_x[dim_x != 1]
  if (length(d) != 2)
    stop("data must be in a matrix form")
  X <- if (length(d) != length(dim_x))
  {
    matrix(X, d[1], d[2])
  }
  else {
    as.matrix(X)
  }
  
  alpha <- 1
  
  X.names = dimnames(X)[[2]]
  if (is.null(X.names))
    X.names = paste("X", 1:ncol(X), sep = "")
  
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names))
    ind.names = 1:nrow(X)
  ## center but scale only if asked
  X <- scale(X, center = TRUE, scale = scale)
  
  svd_mat <- svd(X)
  right_sing_vect <- svd_mat$v
  right_sing_vect <- scale(right_sing_vect, center = TRUE, scale = TRUE)
  n <- nrow(t(X))
  p <- ncol(t(X))
  
  if (ncomp > min(n, p)) {
    message("'ncomp' is too large: reset to ", min(n, p))
    ncomp <- min(n, p)
  }
  if (is.null(w.init))
    w.init <- matrix(1 / sqrt(ncomp), ncomp, ncomp)
  else {
    if (!is.matrix(w.init) || length(w.init) != (ncomp ^ 2))
      stop("w.init is not a matrix or is the wrong size")
  }
  
  X1 <- t(right_sing_vect)[1:ncomp, ]
  
  if (mode == "deflation") {
    unmix_mat <- ica.def(
      X1,
      ncomp,
      tol = tol,
      fun = fun,
      alpha = alpha,
      max.iter = max.iter,
      verbose = FALSE,
      w.init = w.init
    )
  }
  else if (mode == "parallel") {
    unmix_mat <- ica.par(
      X1,
      ncomp,
      tol = tol,
      fun = fun,
      alpha = alpha,
      max.iter = max.iter,
      verbose = FALSE,
      w.init = w.init
    )
  }
  w <- unmix_mat
  independent_mat <- w %*% X1
  #==order independent_mat by kurtosis==#
  kurt <- vector(length = ncomp)
  independent_mat.new <- matrix(nrow = ncomp, ncol = n)
  for (h in 1:ncomp) {
    kurt[h] <-
      (mean(independent_mat[h, ] ^ 4) - 3 * (mean(independent_mat[h, ] ^ 2)) ^
         2)
  }
  for (i in 1:ncomp) {
    independent_mat.new[i, ] <-
      independent_mat[order(kurt, decreasing = TRUE)[i], ]
    independent_mat.new[i, ] <-
      independent_mat.new[i, ] / as.vector(crossprod(independent_mat.new[i, ]))
  }
  
  #== variable selection==#
  v.sparse = matrix(nrow = ncomp, ncol = n)
  for (i in 1:ncomp) {
    nx <- n - keepX[i]
    v.sparse[i,] = ifelse(
      abs(independent_mat.new[i,]) > abs(independent_mat.new[i,][order(abs(independent_mat.new[i,]))][nx]),
      (abs(independent_mat.new[i,]) - abs(independent_mat.new[i,][order(abs(independent_mat.new[i,]))][nx])) * sign(independent_mat.new[i,]),
      0
    )
  }
  independent_mat.new = v.sparse
  
  
  mix_mat <- t(w) %*% solve(w %*% t(w))
  
  ipc_mat = matrix(nrow = p, ncol = ncomp)
  ipc_mat = X %*% t(independent_mat.new)
  ##== force orthogonality ==##
  for (h in 1:ncomp) {
    if (h == 1) {
      ipc_mat[, h] = X %*% (t(independent_mat.new)[, h])
    }
    if (h > 1) {
      ipc_mat[, h] = (lsfit(y = X %*% (t(
        independent_mat.new
      )[, h]), ipc_mat[, 1:(h - 1)], intercept = FALSE)$res)
    }
    ipc_mat[, h] = ipc_mat[, h] / as.vector(sqrt(crossprod(ipc_mat[, h])))
  }
  ##== force over ==##
  # put rownames of loading vectors
  colnames(independent_mat.new) = colnames(X)
  
  result = (
    list(
      X = X,
      ncomp = ncomp,
      keepX = keepX,
      unmixing = t(unmix_mat),
      mixing = t(mix_mat),
      loadings = list(X = t(independent_mat.new)),
      rotation = t(independent_mat.new),
      kurtosis = kurt[order(kurt, decreasing = TRUE)],
      names = list(X = X.names, sample = ind.names)
    )
  )
  
  
  result$x = ipc_mat
  result$variates = list(X = ipc_mat)
  dimnames(result$x) = list(ind.names, paste("IPC", 1:ncol(result$rotation), sep = " "))
  
  class(result) = c("sipca", "ipca", "pca")
  
  #calcul explained variance
  explX = explained_variance(X, result$variates$X, ncomp)
  result$explained_variance = explX
  
  ## keeping this cause we might revert S4 methods
  ## the return us handled in methods now
  ret.call <- FALSE
  if ( isTRUE(ret.call) ) { ## return call
    mcr <- match.call()
    mcr[-1] <- lapply(mcr[-1], eval)
    mcr[[1L]] <- quote(sipca)
    result <- c(call = mcr, result)
  }
  return(invisible(result))
}

## ----------- Generic ----------- 
#' @param ... Aguments passed to the generic.
#' @export
setGeneric("sipca", function(data=NULL, X=NULL, ...) standardGeneric("sipca") )

## ----------- Methods ----------- 

#### Default ####

## in default method data should be NULL, but we make an exception for legacy
## codes where first argument (used to be X) is not named and is taken as data

#' @export
setMethod("sipca", "ANY",function(data=NULL, X=NULL, ...){
  ## evaluate matched call with defaults
  mget(names(formals()), sys.frame(sys.nframe()))
  mc <- match.call()
  mc[-1L] <- lapply(mc[-1L], eval.parent)
  if ( any(class(data) %in% c('matrix', 'data.frame'))) {
    .deprecate_ufma() ## deprecate unnamed first matrix argument
    mc$X <- mc$data
    mc$data <- NULL
    mc[[1L]] <- quote(.sipca)
  } else if ( !(missing(data) || is.null(data)) ) {
    .stop("'data' is not valid, see ?sipca.", .subclass = "inv_data")
  }
  mc$data <- NULL ## not needed for the internal
  mc[[1L]] <- quote(.sipca) ## call to internal
  result <- eval(mc)
  ## include the call if asked
  .call_return(result, mc$ret.call, mcr = match.call(), 
               fun.name = 'sipca')
} )

#### MultiAssayExperiment ####
#' @export
#' @rdname sipca
setMethod("sipca", "MultiAssayExperiment", function(data=NULL, X=NULL, ...){
  ## evaluate matched call with defaults
  mget(names(formals()), sys.frame(sys.nframe()))
  mc <- match.call()
  mc[-1] <- lapply(mc[-1], eval.parent)
  ## pipeline to adjust and send the method args to the internal
  result <- .pcaMethodsHelper(mc, fun = 'sipca')
  ## if asked, append the evaluted args to the output list
  .call_return(result, mc$ret.call,
               mcr = match.call(), fun.name = 'sipca')
})
