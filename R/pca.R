#' Principal Components Analysis
#' 
#' Performs a principal components analysis on the given data matrix that can
#' contain missing values. If data are complete 'pca' uses Singular Value
#' Decomposition, if there are some missing values, it uses the NIPALS
#' algorithm.
#' 
#' The calculation is done either by a singular value decomposition of the
#' (possibly centered and scaled) data matrix, if the data is complete or by
#' using the NIPALS algorithm if there is data missing. Unlike
#' \code{\link{princomp}}, the print method for these objects prints the
#' results in a nice format and the \code{plot} method produces a bar plot of
#' the percentage of variance explaned by the principal components (PCs).
#' 
#' When using NIPALS (missing values), we make the assumption that the first
#' (\code{min(ncol(X),} \code{nrow(X)}) principal components will account for
#' 100 \% of the explained variance.
#' 
#' Note that \code{scale= TRUE} cannot be used if there are zero or constant
#' (for \code{center = TRUE}) variables.
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
#' 
#' @param X a numeric matrix (or data frame) which provides the data for the
#' principal components analysis. It can contain missing values.
#' @param ncomp Integer, if data is complete \code{ncomp} decides the number of
#' components and associated eigenvalues to display from the \code{pcasvd}
#' algorithm and if the data has missing values, \code{ncomp} gives the number
#' of components to keep to perform the reconstitution of the data using the
#' NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X),
#' ncol(X))}
#' @param center (Default=TRUE) Logical, whether the variables should be
#' shifted to be zero centered. Alternatively, a vector of length equal the
#' number of columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param scale (Default=FALSE) Logical indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. The default is
#' \code{FALSE} for consistency with \code{prcomp} function, but in general
#' scaling is advisable. Alternatively, a vector of length equal the number of
#' columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param max.iter Integer, the maximum number of iterations in the NIPALS
#' algorithm.
#' @param tol Positive real, the tolerance used in the NIPALS algorithm.
#' @param logratio (Default='none') one of ('none','CLR','ILR'). Specifies the log ratio
#' transformation to deal with compositional values that may arise from
#' specific normalisation in sequencing data. Default to 'none'
#' @param ilr.offset (Default=0.001) When logratio is set to 'ILR', an offset must be input to
#' avoid infinite value after the logratio transform.
#' @param V Matrix used in the logratio transformation if provided.
#' @param multilevel sample information for multilevel decomposition for
#' repeated measurements.
#' @param reconst (Default=TRUE) Logical. If matrix includes missing values, whether
#' \code{nipals} must perform the reconstitution of the data using the
#' \code{ncomp} components. The components are calculated by imputing
#' the missing values when set to TRUE, otherwise the missing values will be ignored.
#' @return \code{pca} returns a list with class \code{"pca"} and \code{"prcomp"}
#' containing the following components: 
#' \item{call}{The function call.}
#' \item{ncomp}{The number of principal components used.}
#' \item{center}{The centering used.}
#' \item{scale}{The scaling used.}
#' \item{names}{List of row and column names of data.}
#' \item{sdev}{The eigenvalues of the covariance/correlation matrix, though
#' the calculation is actually done with the singular values of the data
#' matrix or by using NIPALS.}
#' \item{rotation}{The matrix of variable loadings (i.e., a matrix whose
#' columns contain the eigenvectors).}
#' \item{x}{The value of the rotated data (the centred (and scaled if
#' requested) data multiplied by the rotation/loadings matrix), also called
#' the principal components.}
#' \item{var.tot}{Total variance in the data.}
#' \item{loadings}{Same as 'rotation' to keep the mixOmics spirit}
#' \item{variates}{Same as 'x' to keep the mixOmics spirit}
#' \item{explained_variance}{Explained variance from the multivariate model,
#'   used for plotIndiv}
#' \item{cum.var}{The cumulative explained variance for components.}
#' \item{X}{The input data matrix.}
#' @author Florian Rohart, Kim-Anh Lê Cao, Ignacio González, Al J Abadi
#' @seealso \code{\link{nipals}}, \code{\link{prcomp}}, \code{\link{biplot}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}} and http://www.mixOmics.org
#' for more details.
#' @references On log ratio transformations: Filzmoser, P., Hron, K., Reimann,
#' C.: Principal component analysis for compositional data with outliers.
#' Environmetrics 20(6), 621-632 (2009) Lê Cao K.-A., Costello ME, Lakis VA,
#' Bartolo, F,Chua XY, Brazeilles R, Rondeau P. MixMC: Multivariate insights
#' into Microbial Communities. PLoS ONE, 11(8): e0160169 (2016). On multilevel
#' decomposition: Westerhuis, J.A., van Velzen, E.J., Hoefsloot, H.C., Smilde,
#' A.K.: Multivariate paired data analysis: multilevel plsda versus oplsda.
#' Metabolomics 6(1), 119-128 (2010) Liquet, B., Lê Cao, K.-A., Hocini, H.,
#' Thiebaut, R.: A novel approach for biomarker selection and the integration
#' of repeated measures experiments from two assays. BMC bioinformatics 13(1),
#' 325 (2012)
#' @keywords algebra
#' @export
#' @example ./examples/pca-examples.R
pca <- function(X,
                ncomp = 2,
                center = TRUE,
                scale = FALSE,
                max.iter = 500,
                tol = 1e-09,
                logratio = c('none','CLR','ILR'),
                ilr.offset = 0.001,
                V = NULL,
                multilevel = NULL,
                reconst = TRUE)
{
    
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                   error = function(e) e)
    
    if ("simpleError" %in% class(err))
        stop(err[[1]], ".", call. = FALSE)
    
    #-- X matrix
    if (is.data.frame(X))
        X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
        stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
        stop("infinite values in 'X'.", call. = FALSE)
    
    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
        X.names = paste("V", 1:ncol(X), sep = "")
    
    ind.names = rownames(X)
    if (is.null(ind.names))
        ind.names = 1:nrow(X)
    
    #-- ncomp
    if (is.null(ncomp))
        ncomp = min(nrow(X), ncol(X))
    
    ncomp = round(ncomp)
    
    if (!is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
        stop("invalid value for 'ncomp'.", call. = FALSE)
    
    if (ncomp > min(ncol(X), nrow(X)))
        stop("use smaller 'ncomp'", call. = FALSE)
    
    #-- log.ratio
    choices = c('CLR', 'ILR', 'none')
    logratio = choices[pmatch(logratio, choices)]
    logratio <- match.arg(logratio)

    if (logratio != "none" && any(X < 0))
        stop("'X' contains negative values, you can not log-transform your data"
        )
    
    #-- cheking center and scale
    if (!is.logical(center))
    {
        if (!is.numeric(center) || (length(center) != ncol(X)))
            stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
                 call. = FALSE)
    }
    
    if (!is.logical(scale))
    {
        if (!is.numeric(scale) || (length(scale) != ncol(X)))
            stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
                 call. = FALSE)
    }
    
    #-- max.iter
    if (is.null(max.iter) || !is.numeric(max.iter) || max.iter < 1 || !is.finite(max.iter))
        stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)
    
    #-- tol
    if (is.null(tol) || !is.numeric(tol) || tol < 0 || !is.finite(tol))
        stop("invalid value for 'tol'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    
    #-----------------------------#
    #-- logratio transformation --#
    
    if (is.null(V) & logratio == "ILR") # back-transformation to clr-space, will be used later to recalculate loadings etc
        V = clr.backtransfo(X)
    
    X = logratio.transfo(X = X, logratio = logratio, offset = if(logratio == "ILR") {ilr.offset} else {0})
    
    #as X may have changed
    if (ncomp > min(ncol(X), nrow(X)))
        stop("use smaller 'ncomp'", call. = FALSE)
    
    #-- logratio transformation --#
    #-----------------------------#
    
    #---------------------------------------------------------------------------#
    #-- multilevel approach ----------------------------------------------------#
    
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
    #-- multilevel approach ----------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    X = scale(X, center = center, scale = scale)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    
    if (any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance.",
             call. = FALSE)
    
    cl = match.call()
    cl[[1]] = as.name('pca')
    result = list(call = cl, X = X, ncomp = ncomp,
                  center = if (is.null(cen)) {FALSE} else {cen},
                  scale = if (is.null(sc)) {FALSE} else {sc},
                  names = list(X = X.names, sample = ind.names))
    

    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    if (any(is.na(X))) {
        res <- nipals(X, ncomp = ncomp, reconst = reconst, max.iter = max.iter, tol = tol)
        res$rotation <- res$p
        result <- c(result, res[c("ncomp", "sdev", "var.tot", "loadings", 
                                  "variates", "rotation", "x", "explained_variance", "cum.var")])
        rm(res)
    } else {
        if (logratio %in% c('CLR', 'none')) {
            #-- if data is complete use singular value decomposition
            #-- borrowed from 'prcomp' function
            res = svd(X, nu = 0)
            result <-
                .add_sdev_rotation_and_comps(
                    result = result,
                    X = X,
                    ncomp = ncomp,
                    sdev = res$d,
                    rotation = res$v,
                    row_names = ind.names,
                    col_names = X.names
                )
        } else {
            # if 'ILR', transform data and then back transform in clr space (from RobCompositions package)
            # data have been transformed above
            res = svd(X, nu = max(1, nrow(X) - 1))
            result$sdev = res$d[1:ncomp] / sqrt(max(1, nrow(X) - 1))  # Note: what differs with RobCompo is that they use: cumsum(eigen(cov(X))$values)/sum(eigen(cov(X))$values)
            # calculate loadings using back transformation to clr-space
            result$rotation = V %*% res$v[, 1:ncomp, drop = FALSE]
            # extract component score from the svd, multiply matrix by vector using diag, NB: this differ from our mixOmics PCA calculations
            # NB: this differ also from Filmoser paper, but ok from their code: scores are unchanged
            result$x = res$u[, 1:ncomp, drop = FALSE] %*% diag(res$d[1:ncomp, drop = FALSE])
            names(result$sdev) = paste("PC", 1:length(result$sdev), sep = "")
            dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
            dimnames(result$x) = list(ind.names, paste("PC", 1:ncol(result$x), sep = ""))
            
        }
        ## add variance stats
        result <- .add_var_stats(result)
        result$X <- X
    }
    
    # output multilevel if needed
    if(!is.null(multilevel))
        result=c(result, list(Xw = Xw, design = multilevel))
    
    class(result) = c("pca","prcomp") 
    if(!is.null(multilevel))
        class(result)=c("mlpca",class(result))
    
    return(invisible(result))
}

#' Add sdev, rotation and components to result
#'
#' @param result Result list. $sdev, $rotation and $x will be added to it.
#' @param X A numeric matrix.
#' @param ncomp Integer, number of components.
#' @param sdev numeric vector - sqrt(eig) from different analyses
#' @param rotation Numeric matrix of rotation from different analyses
#' @param row_names Vector of sample names
#' @param col_names Vector of feature names
#' @noRd
#' @return A modified list
.add_sdev_rotation_and_comps <-
    function(result,
             X,
             ncomp,
             sdev,
             rotation,
             row_names,
             col_names) {
        result$sdev = sdev[seq_len(ncomp)] / sqrt(max(1, nrow(X) - 1))
        names(result$sdev) = paste("PC", seq_along(result$sdev), sep = "")
        result$rotation = rotation[, seq_len(ncomp), drop = FALSE]
        dimnames(result$rotation) = list(col_names, paste("PC", seq_len(ncomp), sep = ""))
        result$x = X %*% result$rotation
        dimnames(result$x) = list(row_names, paste("PC", seq_len(ncol(result$x)), sep = ""))
        result$X <- X
        result
    }


#' Add cumulative and per-component explained variance 
#'
#' @param result A list, output of \code{\link{.add_sdev_rotation_and_comps}}
#'
#' @return A list. Modifies input list by adding total variance, loadings and
#' variates (only for consistency of mixOmics outputs), and per-component and
#' cumulative explained variance.
#' @noRd
.add_var_stats <- function(result) {
    result$var.tot=sum(result$X^2) / max(1, nrow(result$X) - 1)# same as all res$d, or variance after nipals replacement of the missing values
    
    # to be similar to other methods, add loadings and variates as outputs
    result$loadings = list(X=result$rotation)
    result$variates = list(X=result$x)
    
    # calculate explained variance
    explained_variance <- with(result, sdev^2 / var.tot)
    result$cum.var = cumsum(explained_variance)
    result$explained_variance = list(X = explained_variance)
    result$X <- NULL
    result
}

