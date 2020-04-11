# =============================================================================
# Calculate the explained variance of one dataset based on its variates
# =============================================================================

#' Calculation of explained variance
#' 
#' This function calculates the variance explained by variates.
#' 
#' 
#' \code{explained_variance} calculates the explained variance of each variates
#' out of the total variance in \code{data}.
#' 
#' @param data numeric matrix of predictors
#' @param variates variates as obtained from a \code{pls} object for instance
#' @param ncomp number of components. Should be lower than the number of
#' columns of \code{variates}
#' @return \code{explained_variance} simply returns the explained variance for
#' each variate.
#' @author Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{spls}}, \code{\link{splsda}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{cim}}, \code{\link{network}}.
#' @keywords regression multivariate
#' @examples
#' 
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' 
#' toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))
#' 
#' ex = explained_variance(toxicity.spls$X, toxicity.spls$variates$X, ncomp =2)
#' 
#' # ex should be the same as
#' toxicity.spls$explained_variance$X
#' 
#' @export
explained_variance = function(data, variates, ncomp)
{
    #check input data
    check = Check.entry.single(data, ncomp)
    data = check$X
    ncomp = check$ncomp
    
    if (anyNA(data))
    {
        warning("NA values put to zero, results will differ from PCA methods
        used with NIPALS")
        isna = is.na(data)
        data[isna] = 0
    }
    nor2x <- sum((data)^2) # total variance in the data
    
	exp.varX = NULL
	for (h in 1:ncomp)
	{
        a <- t(variates[, h, drop=FALSE]) %*% data
        ta = t(a)
        exp_var_new <- a%*%ta /crossprod(variates[, h],variates[, h])/nor2x


	    exp.varX = append(exp.varX, exp_var_new)
	
	}
    names(exp.varX) = paste("comp", 1:ncomp)
    
    # result: vector of length ncomp with the explained variance per component
    exp.varX
}

