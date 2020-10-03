#' Non-linear Iterative Partial Least Squares (NIPALS) algorithm
#' 
#' This function performs NIPALS algorithm, i.e. the singular-value
#' decomposition (SVD) of a data table that can contain missing values.
#' 
#' The NIPALS algorithm (Non-linear Iterative Partial Least Squares) has been
#' developed by H. Wold at first for PCA and later-on for PLS. It is the most
#' commonly used method for calculating the principal components of a data set.
#' It gives more numerically accurate results when compared with the SVD of the
#' covariance matrix, but is slower to calculate.
#' 
#' This algorithm allows to realize SVD with missing data, without having to
#' delete the rows with missing data or to estimate the missing data.
#' 
#' @inheritParams pca
#' @param reconst logical that specify if \code{nipals} must perform the
#' reconstitution of the data using the \code{ncomp} components.
#' @return An object of class 'mixo_nipals' contaning slots: 
#' \item{call}{The function call.}
#' \item{eig}{Vector containing the pseudo-singular values of \code{X}, of length
#' \code{ncomp}.}
#' \item{p}{Matrix whose columns contain the right singular vectors of \code{X}.}
#' \item{t}{Matrix whose columns contain the left singular vectors of \code{X}.
#' Note that for a complete data matrix X, the return values \code{eig},
#' \code{t} and \code{p} such that \code{X = t * diag(eig) * t(p)}.}
#' \item{ncomp}{The number of principal components used.}
#' \item{rec}{If \code{reonst=TRUE}, matrix obtained by the reconstitution of
#' the data using the \code{ncomp} components.}
#' \item{sdev}{Same as 'eig' - for mixOmics consistency.}
#' \item{var.tot}{Total variance in the data.}
#' \item{loadings}{ame as 'p' to keep the mixOmics spirit}
#' \item{x}{the value of the rotated data (the centred (and scaled if
#' requested) data multiplied by the rotation/loadings matrix), also called
#' the principal components.}
#' \item{variates}{Same as 'x' to keep the mixOmics spirit}
#' \item{explained_variance}{explained variance of each component.}
#' \item{cum.var}{The cumulative explained variance for components.}
#'  
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Le Cao, Al J Abadi
#' @seealso \code{\link{svd}}, \code{\link{princomp}}, \code{\link{prcomp}},
#' \code{\link{eigen}} and http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' Wold H. (1975). Path models with latent variables: The NIPALS approach. In:
#' Blalock H. M. et al. (editors). \emph{Quantitative Sociology: International
#' perspectives on mathematical and statistical model building}. Academic
#' Press, N.Y., 307-357.
#' @keywords algebra multivariate
#' @export
#' @examples
#' 
#' ## Hilbert matrix
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#' X.na <- X <- hilbert(9)[, 1:6]
#' 
#' ## Hilbert matrix with missing data
#' idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
#' X.na[idx.na == 0] <- NA
#' X.rec <- nipals(X.na, reconst = TRUE)$rec
#' round(X, 2)
#' round(X.rec, 2)
nipals <- function (X,
                    ncomp = 2,
                    reconst = FALSE,
                    max.iter = 500,
                    tol = 1e-06)
{
    #-- X matrix
    if (is.data.frame(X))
        X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
        stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
        stop("infinite values in 'X'.", call. = FALSE)
    
    nc = ncol(X)
    nr = nrow(X)
    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
        X.names = paste("V", 1:ncol(X), sep = "")
    
    ind.names = rownames(X)
    if (is.null(ind.names))
        ind.names = 1:nrow(X)
    
    #-- ncomp
    if (is.null(ncomp) ||
        !is.numeric(ncomp) ||
        ncomp < 1 || !is.finite(ncomp) || ncomp > min(nr, nc))
        stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    #-- reconst
    if (!is.logical(reconst))
        stop("'reconst' must be a logical constant (TRUE or FALSE).",
             call. = FALSE)
    
    #-- max.iter
    if (is.null(max.iter) || max.iter < 1 || !is.finite(max.iter))
        stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)
    
    #-- tol
    if (is.null(tol) || tol < 0 || !is.finite(tol))
        stop("invalid value for 'tol'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    #-- initialisation des matrices --#
    X.input <- X
    p = matrix(nrow = nc, ncol = ncomp)
    t.mat = matrix(nrow = nr, ncol = ncomp)
    eig = vector("numeric", length = ncomp)
    nc.ones = rep(1, nc)
    nr.ones = rep(1, nr)
    is.na.X = is.na(X.input)
    na.X <- ifelse(any(is.na.X), TRUE, FALSE)
    #-- boucle sur h --#
    for (h in 1:ncomp)
    {
        th = X[, which.max(apply(X, 2, var, na.rm = TRUE))]
        if (any(is.na(th))) th[is.na(th)] = 0
        ph.old = rep(1 / sqrt(nc), nc)
        ph.new = vector("numeric", length = nc)
        iter = 1
        diff = 1
        
        if (na.X)
        {
            X.aux = X
            X.aux[is.na.X] = 0
        }
        
        while (diff > tol & iter <= max.iter)
        {
            if (na.X)
            {
                ph.new = crossprod(X.aux, th)
                Th = drop(th) %o% nc.ones
                Th[is.na.X] = 0
                th.cross = crossprod(Th)
                ph.new = ph.new / diag(th.cross)
            } else {
                ph.new = crossprod(X, th) / drop(crossprod(th))
            }
            
            ph.new = ph.new / drop(sqrt(crossprod(ph.new)))
            
            if (na.X)
            {
                th = X.aux %*% ph.new
                P = drop(ph.new) %o% nr.ones
                P[t(is.na.X)] = 0
                ph.cross = crossprod(P)
                th = th / diag(ph.cross)
            } else {
                th = X %*% ph.new / drop(crossprod(ph.new))
            }
            
            diff = drop(sum((ph.new - ph.old)^2, na.rm = TRUE))
            ph.old = ph.new
            iter = iter + 1
        }
        
        if (iter > max.iter)
            warning(paste("Maximum number of iterations reached for comp.", h))
        
        X = X - th %*% t(ph.new)
        p[, h] = ph.new
        t.mat[, h] = th
        eig[h] = sum(th * th, na.rm = TRUE)
    }
    
    eig = sqrt(eig)
    t.mat = scale(t.mat, center = FALSE, scale = eig)
    attr(t.mat, "scaled:scale") = NULL
    result = list(call = match.call(), eig = eig, p = p, t = t.mat, ncomp = ncomp)
    
    if (reconst)
    {
        if (ncomp < 5)
            message("\nconsider high 'ncomp' for more accurate ",
                    "imputation of the missing values.\n")
        X.hat = t.mat %*% diag(eig) %*% t(p)
        
        colnames(X.hat) = colnames(X)
        rownames(X.hat) = rownames(X)
        result$rec <- X.hat
    }
    
    if (isFALSE(reconst)) { ## replace NA with 0 with messages
        message('\nreplacing missing values with 0 for calculation',
            'of components as `reconst=FALSE` used.')
        if (is.null(attr(X.input, 'scaled:center'))) {
            message('This will have implications particularly if data are not centred.',
                    ' (we recommend the use of `center = TRUE`).\n')
        }
        X.input[is.na.X] <- 0
    } else {
        X.input[is.na.X] <- result$rec[is.na.X]
    }
    result <-
        .add_sdev_rotation_and_comps(
            result = result,
            X = X.input,
            ncomp = ncomp,
            sdev = result$eig,
            rotation = result$p,
            row_names = ind.names,
            col_names = X.names
        )
    ## add variance stats
    result <- .add_var_stats(result)
    class(result) <- c('pca', 'mixo_nipals')
    return(invisible(result))
}
