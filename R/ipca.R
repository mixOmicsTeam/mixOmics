####################################################################################
## ---------- internal
.ipca = function(X,
                 ncomp = 2,
                 mode = c("deflation", "parallel"),
                 fun = c("logcosh", "exp", "kur"),
                 scale = FALSE,
                 w.init = NULL,
                 max.iter = 200,
                 tol = 1e-04) {

  ## ----------------------------------- checks
  mode <- .matchArg(mode)
  fun <- .matchArg(fun)

  mcd <- mget(names(formals()),sys.frame(sys.nframe())) ## match.call and defaults
  err = tryCatch(mcd, error = function(e) e) ## see if arguments can be evaluated
  if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)

  ## check pca entries for match.call and defaults and do necessary adjustments to
  ## arguments in this environement from within function
  .pcaEntryChecker(mcd, check.keepX = FALSE, check.center = FALSE, check.NA = TRUE)

  nc = ncol(X)
  nr = nrow(X)

  #-- put a names on the rows and columns of X --#
  X.names = colnames(X)
  if (is.null(X.names))
    X.names = paste("V", 1:nc, sep = "")

  ind.names = rownames(X)
  if (is.null(ind.names))
    ind.names = 1:nr

## ----------------------------------- scale
  X = scale(X, center = TRUE, scale = scale)
  sc = attr(X, "scaled:scale")

  if (any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance.",
         call. = FALSE)

  #-- mode
  choices = c("deflation", "parallel")
  mode = choices[pmatch(mode, choices)]

  if (is.na(mode))
    stop("'mode' should be one of 'deflation' or 'parallel'.", call. = FALSE)

  #-- fun
  choices = c("logcosh", "exp")
  fun = choices[pmatch(fun, choices)]

  if (is.na(fun))
    stop("'fun' should be one of 'logcosh' or 'exp'.", call. = FALSE)

  #-- w.init
  if (is.null(w.init))
  {
    w.init = matrix(1 / sqrt(ncomp), ncomp, ncomp)
  } else {
    if (!is.matrix(w.init) ||
        length(w.init) != (ncomp ^ 2) || !is.numeric(w.init))
      stop("'w.init' is not a numeric matrix or is the wrong size", call. = FALSE)
  }

  if (any(is.infinite(w.init)))
    stop("infinite values in 'w.init'.", call. = FALSE)

  if (sum(w.init == 0, na.rm = TRUE) == length(w.init))
    stop("'w.init' has to be a non-zero matrix", call. = FALSE)

  if (any(is.na(w.init)))
    stop("'w.init' has to be a numeric matrix matrix with non-NA values",
         call. = FALSE)


  #-- max.iter
  if (is.null(max.iter) ||
      !is.numeric(max.iter) || max.iter < 1 || !is.finite(max.iter))
    stop("invalid value for 'max.iter'.", call. = FALSE)

  max.iter = round(max.iter)

  #-- tol
  if (is.null(tol) ||
      !is.numeric(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)


  #-- end checking --#
  #------------------#

  #-- ipca approach ----------------------------------------------------------#
  #---------------------------------------------------------------------------#

  V = svd(X)$v
  V = scale(V, center = TRUE, scale = TRUE)
  V = t(V)[1:ncomp, , drop = FALSE]

  if (mode == "deflation")
  {
    W = ica.def(
      V,
      ncomp = ncomp,
      tol = tol,
      fun = fun,
      alpha = 1,
      max.iter = max.iter,
      verbose = FALSE,
      w.init = w.init
    )
  } else if (mode == "parallel") {
    W = ica.par(
      V,
      ncomp = ncomp,
      tol = tol,
      fun = fun,
      alpha = 1,
      max.iter = max.iter,
      verbose = FALSE,
      w.init = w.init
    )
  }

  #-- independent loadings --#
  S = matrix(W %*% V, nrow = ncomp)

  #-- order independent loadings by kurtosis --#
  kurt = apply(S, 1, function(x) {
    n = length(x)
    x = x - mean(x)
    n * sum(x ^ 4) / (sum(x ^ 2) ^ 2) - 3
  })
  ord = order(kurt, decreasing = TRUE)
  kurt = kurt[ord]
  S = S[ord, , drop = FALSE]
  norm = apply(S, 1, function(x) {
    crossprod(x)
  })
  S = t(sweep(S, 1, norm, "/"))

  #-- independent PCs / force orthonormality --#
  ipc = matrix(nrow = nr, ncol = ncomp)
  ipc[, 1] = X %*% S[, 1]
  ipc[, 1] = ipc[, 1] / as.numeric(sqrt(crossprod(ipc[, 1])))

  if (ncomp > 1)
  {
    for (h in 2:ncomp)
    {
      ipc[, h] = lsfit(y = X %*% S[, h], ipc[, 1:(h - 1)], intercept = FALSE)$res
      ipc[, h] = ipc[, h] / as.numeric(sqrt(crossprod(ipc[, h])))
    }
  }


  #-- output -----------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  dimnames(S) = list(X.names, paste("IPC", 1:ncol(S), sep = ""))
  dimnames(ipc) = list(ind.names, paste("IPC", 1:ncol(ipc), sep = ""))

  {
    ## keeping this for the exported generic's benefit, but is not necessary as far as methods are concerned
    mcr = match.call() ## match call for returning
    mcr[[1]] = as.name('ipca')
    mcr[-1L] <- lapply(mcr[-1L], eval.parent)
  }

  result = list(
    call = mcr,
    X = X,
    ncomp = ncomp,
    x = ipc,
    loadings = list(X = S),
    rotation = S,
    variates = list(X = ipc),
    kurtosis = kurt,
    unmixing = t(W),
    mixing = t(t(W) %*% solve(W %*% t(W))),
    names = list(var = X.names, sample = ind.names)
  )

  class(result) = c("ipca", "pca")

  #calcul explained variance
  explX = explained_variance(X, result$variates$X, ncomp)
  result$explained_variance = explX

  return(invisible(result))
}


#' @title Independent Principal Component Analysis
#'
#' @description  Performs independent principal component analysis on the given data matrix,
#' a combination of Principal Component Analysis and Independent Component
#' Analysis.
#'
#' In PCA, the loading vectors indicate the importance of the variables in the
#' principal components. In large biological data sets, the loading vectors
#' should only assign large weights to important variables (genes, metabolites
#' ...). That means the distribution of any loading vector should be
#' super-Gaussian: most of the weights are very close to zero while only a few
#' have large (absolute) values.
#'
#' However, due to the existence of noise, the distribution of any loading
#' vector is distorted and tends toward a Gaussian distribtion according to the
#' Central Limit Theroem. By maximizing the non-Gaussianity of the loading
#' vectors using FastICA, we obtain more noiseless loading vectors. We then
#' project the original data matrix on these noiseless loading vectors, to
#' obtain independent principal components, which should be also more noiseless
#' and be able to better cluster the samples according to the biological
#' treatment (note, IPCA is an unsupervised approach).
#'
#' \bold{Algorithm} 1. The original data matrix is centered.
#'
#' 2. PCA is used to reduce dimension and generate the loading vectors.
#'
#' 3. ICA (FastICA) is implemented on the loading vectors to generate
#' independent loading vectors.
#'
#' 4. The centered data matrix is projected on the independent loading vectors
#' to obtain the independent principal components.
#'
## ----------------------------------- Parameters
#' @inheritParams pca
#' @param mode character string. What type of algorithm to use when estimating
#' the unmixing matrix, choose one of \code{"deflation"}, \code{"parallel"}.
#' @param fun the function used in approximation to neg-entropy in the FastICA
#' algorithm. Default set to \code{logcosh}, see details of FastICA.
#' @param w.init initial un-mixing matrix (unlike FastICA, this matrix is fixed
#' here).

## ----------------------------------- Value
#' @return \code{ipca} returns a list with class \code{"ipca"} containing the
#' following components: \item{ncomp}{the number of independent principal
#' components used.} \item{unmixing}{the unmixing matrix of size (ncomp x
#' ncomp)} \item{mixing}{the mixing matrix of size (ncomp x ncomp)}
#' \item{X}{the centered data matrix} \item{x}{the indepenent principal
#' components} \item{loadings}{the independent loading vectors}
#' \item{kurtosis}{the kurtosis measure of the independent loading vectors}

## ----------------------------------- Misc
#' @author Fangzhou Yao and Jeff Coquery.
#' @seealso \code{\link{sipca}}, \code{\link{pca}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, and http://www.mixOmics.org for more details.
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
## ----------------------------------- Examples
#' @example examples/ipca-example.R

####################################################################################
## ---------- Generic
#' @param ... Aguments passed to the generic.
#' @export
ipca <- function(X=NULL, data=NULL, ncomp=2, keepX=NULL, ...) UseMethod('ipca')

####################################################################################
## ---------- Methods

## ----------------------------------- X=matrix
#' @rdname ipca
#' @export
ipca.default <- .ipca

## ----------------------------------- X= assay name from data
#' @importFrom SummarizedExperiment assay
#' @rdname ipca
#' @export
ipca.character <- function(X=NULL, data=NULL, ncomp=2, keepX=NULL, ...){
  .pcaMethodsHelper(match.call(), fun = 'ipca')
}
