# =============================================================================
# get.confusion_matrix: create confusion table between a vector of true classes
#   and a vector of predicted classes
# =============================================================================

# truth: factor of true classes
# levels: levels of the 'truth' factor. Optional parameters that can be used
#   when there are some missing levels in `truth' compared to the fitted model
# predicted: vector of predicted classes. Can contain NA.

#' Create confusion table and calculate the Balanced Error Rate
#' 
#' Create confusion table between a vector of true classes and a vector of
#' predicted classes, calculate the Balanced Error rate
#' 
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#' 
#' get.BER(confusion)
#' @param truth A factor vector indicating the true classes of the samples
#' (typically \code{Y} from the training set).
#' @param all.levels Levels of the 'truth' factor. Optional parameter if there
#' are some missing levels in \code{truth} compared to the fitted predicted
#' model
#' @param predicted Vector of predicted classes (typically the prediction from
#' the test set). Can contain NA.
#' @param confusion result from a \code{get.confusion_matrix} to calculate the
#' Balanced Error Rate
#' @return \code{get.confusion_matrix} returns a confusion matrix.
#' 
#' \code{get.BER} returns the BER from a confusion matrix
#' @author Florian Rohart, Al J Abadi
#' @seealso \code{\link{predict}}.
#' @references
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @examples
#' 
#' 
#' # Example
#' # -----------------------------------
#' 
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- as.factor(liver.toxicity$treatment[, 4])
#' 
#' ## if training is perfomed on 4/5th of the original data
#' samp <- sample(1:5, nrow(X), replace = TRUE)
#' test <- which(samp == 1)   # testing on the first fold
#' train <- setdiff(1:nrow(X), test)
#' 
#' plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
#' test.predict <- predict(plsda.train, X[test, ], dist = "max.dist")
#' Prediction <- test.predict$class$max.dist[, 2]
#' 
#' # the confusion table compares the real subtypes with the predicted subtypes for a 2 component model
#' confusion.mat = get.confusion_matrix(truth = Y[test],
#' predicted = Prediction)
#' 
#' get.BER(confusion.mat)
#' 
#' @export
get.confusion_matrix = function(truth, all.levels, predicted)
{
    
    
    if (length(truth) != length(predicted))
        stop("'truth' and 'predicted' must be of same length")
    
    if (!is.factor(truth))
        truth = factor(truth)
    
    if (missing(all.levels))
        all.levels = levels(truth)
    
    #print(all.levels)
    
    nlevels.truth = length(all.levels)
    
    ClassifResult = array(0, c(nlevels.truth, nlevels.truth + 1))
    #+1 for NA prediction
    rownames(ClassifResult) = all.levels
    colnames(ClassifResult) = paste("predicted.as.", c(all.levels, "NA"),
                                    sep = "")
    
    #--------record of the classification accuracy for each level of Y
    for (i in 1:nlevels.truth)
    {
        ind.i = which(truth == all.levels[i])
        for (ij in 1:nlevels.truth)
            ClassifResult[i, ij] = sum(predicted[ind.i] == all.levels[ij],
                                       na.rm = TRUE)
        
        # if some NA, we add them in the last column (ij+1 = nlevels.truth + 1)
        if (sum(is.na(predicted[ind.i])) > 0)
            ClassifResult[i, ij + 1] =  sum(is.na(predicted[ind.i]))
    }
    
    # if no NA in the prediction, we remove the last column
    if (sum(is.na(predicted)) == 0)
        ClassifResult = ClassifResult [, -(nlevels.truth + 1)]
    
    ClassifResult
}

#' @noRd
#' @export
# calculate BER from a confusion matrix
get.BER = function(confusion)
{
    
    #if(!is.numeric(X)| !is.matrix(X) | length(dim(X)) != 2 | nrow(X)!=ncol(X))
    #stop("'X' must be a square numeric matrix")
    
    nlev = nrow(confusion)
    #calculation of the BER
    ClassifResult.temp = confusion
    diag(ClassifResult.temp) = 0
    BER = sum(
        apply(ClassifResult.temp, 1, sum, na.rm = TRUE) / apply(confusion, 1,
                                                                sum, na.rm = TRUE),
        na.rm = TRUE
    ) / nlev
    return(BER)
}

