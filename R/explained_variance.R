# =============================================================================
# Calculate variance explained in one dataset based on its components
# =============================================================================

#' Calculation of explained variance
#' 
#' This function calculates the variance explained by their own variates (components) based on redundancy.
#' 
#' 
#' \code{explained_variance} calculates the variance explained by each variate / component and 
#' divides by the total variance in \code{data} using the definition of 'redundancy'. This applies to 
#' any component-based approaches.
#' 
#' @param data numeric matrix of predictors
#' @param variates variates as obtained from a \code{pls} object for instance
#' @param ncomp number of components. Should be lower than the number of
#' columns of \code{variates}
#' @return \code{explained_variance} returns the explained variance for
#' each variate.
#' @details Variance explained by component \eqn{t_h} in \eqn{X} for dimension \eqn{h}:
#' \deqn{Rd(X, t_h) = \frac{1}{p} \sum_{j = 1}^p \mbox{cor}^2(X^j, t_h)}
#' where \eqn{X^j} is the variable centered and scaled, \eqn{p} is the total number of variables.
#' @references
#' Tenenhaus, M.,  La Régression PLS théorie et pratique (1998). Technip, Paris, chap2.
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
explained_variance <- function(data, variates, ncomp)
{
  #check input data
  check = Check.entry.single(data, ncomp)
  data = check$X
  ncomp = check$ncomp
  
  if (anyNA(data))
  {
    warning("NA values are set to zero to estimate the amount of explained variance")
    isna = is.na(data)
    data[isna] = 0
  }
  nor2x <- sum((data)^2) # total variance in the data
  
  exp.varX = NULL
  for (h in 1:ncomp)
  {
    a <- t(variates[, h, drop=FALSE]) %*% data
    ta = t(a)
    # this is equivalent to calculate the redundancy as detailed in the help file
    exp_var_new <- a%*%ta /crossprod(variates[, h],variates[, h])/nor2x
    
    if (anyNA(data))
    {
        warning("\nNA values set to zero for explained variance calculations")
        isna = is.na(data)
        data[isna] = 0
    }
    nor2x <- sum((data)^2) # total variance in the data
    
    exp.varX = append(exp.varX, exp_var_new)
    
  }
  names(exp.varX) = paste("comp", 1:ncomp)
  
  # result: vector of length ncomp with the explained variance per component
  exp.varX
}
