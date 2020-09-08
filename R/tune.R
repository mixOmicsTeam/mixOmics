# ========================================================================================================
# tun: chose the optimal number of parameters per component on a "method"
# ========================================================================================================

#' Generic function to choose the parameters in the different methods in
#' mixOmics
#' 
#' Wrapper of all tuning functions.
#' 
#' The \code{tune} function called the function \code{predict}. more details
#' about most arguments are detailed in \code{?predict}.
#' 
#' Also see the help file corresponding to your \code{method}, e.g.
#' \code{tune.splsda}. Note that only the arguments used in the tune function
#' corresponding to \code{method} are passed on.
#' 
#' Some details on the use of the nrepeat argument are provided in
#' \code{?perf}.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017). More
#' details about the PLS modes are in \code{?pls}.
#' 
#' @param method This parameter is used to pass all other argument to the
#' suitable function. \code{method} has to be one of the following: "spls",
#' "splsda", "mint.splsda", "rcc", "pca".
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y Either a factor or a class vector for the discrete outcome, or a
#' numeric vector or matrix of continuous responses (for multi-response
#' models).
#' @param multilevel Design matrix for multilevel anaylis (for repeated
#' measurements) that indicates the repeated measures on each individual, i.e.
#' the individuals ID. See Details.
#' @param ncomp the number of components to include in the model.
#' @param study grouping factor indicating which samples are from the same
#' study
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set
#' @param test.keepY If \code{method = 'spls'}, numeric vector for the
#' different number of variables to test from the \eqn{Y} data set
#' @param already.tested.X Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{X} data set on
#' the firsts components.
#' @param already.tested.Y if \code{method = 'spls'} and \code{if(ncomp > 1)}
#' numeric vector indicating the number of variables to select from the \eqn{Y}
#' data set on the first components
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' @param grid1,grid2 vector numeric defining the values of \code{lambda1} and
#' \code{lambda2} at which cross-validation score should be computed. Defaults
#' to \code{grid1=grid2=seq(0.001, 1, length=5)}.
#' @param validation character.  What kind of (internal) validation to use,
#' matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is
#' \code{"Mfold"}.
#' @param folds the folds in the Mfold cross-validation. See Details.
#' @param dist distance metric to estimate the
#' classification error rate, should be a subset of \code{"centroids.dist"},
#' \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).
#' @param measure Two misclassification measure are available: overall
#' misclassification error \code{overall} or the Balanced Error Rate \code{BER}
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model.
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default value is FALSE
#' @param logratio one of ('none','CLR'). Default to 'none'
#' @param center a logical value indicating whether the variables should be
#' shifted to be zero centered. Alternately, a vector of length equal the
#' number of columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. The default is
#' \code{FALSE} for consistency with \code{prcomp} function, but in general
#' scaling is advisable. Alternatively, a vector of length equal the number of
#' columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param max.iter Integer, the maximum number of iterations.
#' @param tol Numeric, convergence tolerance criteria.
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @param cpus Integer, number of cores to use for parallel processing. 
#' Currently only available for \code{method = "spls"}
#' @return Depending on the type of analysis performed and the input arguments,
#' a list that may contain:
#' 
#' \item{error.rate}{returns the prediction error for each \code{test.keepX} on
#' each component, averaged across all repeats and subsampling folds. Standard
#' deviation is also output. All error rates are also available as a list.}
#' \item{choice.keepX}{returns the number of variables selected (optimal keepX)
#' on each component.} \item{choice.ncomp}{For supervised models; returns the
#' optimal number of components for the model for each prediction distance
#' using one-sided t-tests that test for a significant difference in the mean
#' error rate (gain in prediction) when components are added to the model. See
#' more details in Rohart et al 2017 Suppl. For more than one block, an optimal
#' ncomp is returned for each prediction framework.}
#' \item{error.rate.class}{returns the error rate for each level of \code{Y}
#' and for each component computed with the optimal keepX}
#' 
#' \item{predict}{Prediction values for each sample, each \code{test.keepX},
#' each comp and each repeat. Only if light.output=FALSE}
#' \item{class}{Predicted class for each sample, each \code{test.keepX}, each
#' comp and each repeat. Only if light.output=FALSE}
#' 
#' \item{auc}{AUC mean and standard deviation if the number of categories in
#' \code{Y} is greater than 2, see details above. Only if auc = TRUE}
#' 
#' \item{cor.value}{only if multilevel analysis with 2 factors: correlation
#' between latent variables.}
#' @author Florian Rohart, Francois Bartolo, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{tune.rcc}}, \code{\link{tune.mint.splsda}},
#' \code{\link{tune.pca}}, \code{\link{tune.splsda}},
#' \code{\link{tune.splslevel}} and http://www.mixOmics.org for more details.
#' @references 
#' Singh A., Shannon C., Gautier B., Rohart F., Vacher M., Tebbutt S.
#' and Lê Cao K.A. (2019), DIABLO: an integrative approach for identifying key 
#' molecular drivers from multi-omics assays, Bioinformatics, 
#' Volume 35, Issue 17, 1 September 2019, Pages 3055–3062.	
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' MINT:
#' 
#' Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017). MINT: A
#' multivariate integrative approach to identify a reproducible biomarker
#' signature across multiple experiments and platforms. BMC Bioinformatics
#' 18:128.
#' 
#' PLS and PLS citeria for PLS regression: Tenenhaus, M. (1998). \emph{La
#' regression PLS: theorie et pratique}. Paris: Editions Technic.
#' 
#' Chavent, Marie and Patouille, Brigitte (2003). Calcul des coefficients de
#' regression et du PRESS en regression PLS1. \emph{Modulad n}, \bold{30} 1-11.
#' (this is the formula we use to calculate the Q2 in perf.pls and perf.spls)
#' 
#' Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction
#' (MSEP) Estimates for Principal Component Regression (PCR) and Partial Least
#' Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9),
#' 422-429.
#' 
#' sparse PLS regression mode:
#' 
#' Lê Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#' 
#' One-sided t-tests (suppl material):
#' 
#' Rohart F, Mason EA, Matigian N, Mosbergen R, Korn O, Chen T, Butcher S,
#' Patel J, Atkinson K, Khosrotehrani K, Fisk NM, Lê Cao K-A&, Wells CA&
#' (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ
#' 4:e1845.
#' @keywords regression multivariate
#' @export
#' @example ./examples/tune-examples.R
tune <-
    function (method,
              # choice of "spls", "splsda", "mint.splsda", "rcc", "pca"
              X,
              Y,
              multilevel = NULL,
              ncomp,
              study,
              # mint.splsda
              test.keepX = c(5, 10, 15),
              # all but pca, rcc
              test.keepY = NULL,
              # rcc, multilevel
              already.tested.X,
              # all but pca, rcc
              already.tested.Y,
              #multilevel
              mode = c("regression", "canonical", "invariant", "classic"),
              # multilevel
              nrepeat = 1,
              #multilevel, splsda
              grid1 = seq(0.001, 1, length = 5),
              # rcc
              grid2 = seq(0.001, 1, length = 5),
              # rcc
              validation = "Mfold",
              # all but pca
              folds = 10,
              # all but pca
              dist = "max.dist",
              # all but pca, rcc
              measure = ifelse(method == "spls", "MSE", "BER"),
              # all but pca, rcc
              auc = FALSE,
              progressBar = FALSE,
              # all but pca, rcc
              near.zero.var = FALSE,
              # all but pca, rcc
              logratio = c('none','CLR'),
              # all but pca, rcc
              center = TRUE,
              # pca
              scale = TRUE,
              # mint, splsda
              max.iter = 100,
              #pca
              tol = 1e-09,
              #pca
              light.output = TRUE,
              # mint, splsda
              cpus = 1
              
    )
    {
        choice.method = c("spls", "splsda", "mint.splsda", "rcc", "pca")
        method = match.arg(method, choice.method)
        mode <- match.arg(mode)
        
        if (method == "mint.splsda") {
            message("Calling 'tune.mint.splsda' with Leave-One-Group-Out Cross Validation (nrepeat = 1)")
            
            if (missing(ncomp))
                ncomp = 1
            
            result = tune.mint.splsda(X = X, Y = Y,
                                      ncomp = ncomp,
                                      study = study,
                                      test.keepX = test.keepX,
                                      already.tested.X = already.tested.X,
                                      dist = dist,
                                      measure = measure,
                                      auc = auc,
                                      progressBar = progressBar,
                                      scale = scale,
                                      tol = tol,
                                      max.iter = max.iter,
                                      near.zero.var = near.zero.var,
                                      light.output = light.output)
            
        } else if (method == "rcc") {
            message("Calling 'tune.rcc'")
            
            result = tune.rcc(X = X,
                              Y = Y,
                              grid1 = grid1,
                              grid2 = grid2,
                              validation = validation,
                              folds = folds,
                              plot = plot)
            
        } else if (method == "pca") {
            message("Calling 'tune.pca'")
            
            if (missing(ncomp))
                ncomp = NULL
            
            result = tune.pca(X = X,
                              ncomp = ncomp,
                              center = center,
                              scale = scale,
                              max.iter = max.iter,
                              tol = tol)
            
            
        } else if (method == "splsda") {
            
            message("Calling 'tune.splsda'")
            
            if (missing(ncomp))
                ncomp = 1
            
            result = tune.splsda (X = X, Y = Y,
                                  ncomp = ncomp,
                                  test.keepX = test.keepX,
                                  already.tested.X = already.tested.X,
                                  validation = validation,
                                  folds = folds,
                                  dist = dist ,
                                  measure = measure,
                                  auc = auc,
                                  progressBar = progressBar,
                                  max.iter = max.iter,
                                  near.zero.var = near.zero.var,
                                  nrepeat = nrepeat,
                                  logratio = logratio,
                                  multilevel = multilevel,
                                  light.output = light.output,
                                  cpus = cpus)
        } else if (method == "spls") {
            if(missing(multilevel))
            {
                message("Calling 'tune.spls'")
                
                result = tune.spls(X = X, Y = Y, ncomp = ncomp, test.keepX = test.keepX, 
                                   already.tested.X = already.tested.X, validation = validation, 
                                   folds = folds, measure = measure,scale = scale, 
                                   progressBar = progressBar, tol = tol, max.iter = max.iter, 
                                   near.zero.var = near.zero.var, nrepeat = nrepeat, 
                                   light.output = light.output, cpus = cpus)
            } else {
                message("Calling 'tune.splslevel' with method = 'spls'")
                
                if (missing(ncomp))
                    ncomp = 1
                if (missing(already.tested.Y))
                    already.tested.Y = NULL
                
                result = tune.splslevel(X = X, Y = Y,
                                        multilevel = multilevel,
                                        mode = mode,
                                        ncomp = ncomp, test.keepX = test.keepX, test.keepY = test.keepY,
                                        already.tested.X = already.tested.X, already.tested.Y = already.tested.Y)
            }
        }
        
        result$call = match.call()
        return(result)
    }
