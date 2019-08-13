#############################################################################################################
# Authors:
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 25-05-2016
# last modified: 24-08-2016
#
# Copyright (C) 2016
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################

#---------------------------------------------------
# perf for mint.spls and mint.pls object
#---------------------------------------------------

#' Predict Method for (mint).(block).(s)pls(da) methods
#'
#' Predicted values based on PLS models. New responses and variates are
#' predicted using a fitted model and a new matrix of observations.
#'
#' \code{predict} produces predicted values, obtained by evaluating the
#' PLS-derived methods, returned by \code{(mint).(block).(s)pls(da)} in the
#' frame \code{newdata}. Variates for \code{newdata} are also returned. Please
#' note that this method performs multilevel decomposition and/or log ratio
#' transformations if needed (\code{multilevel} is an input parameter while
#' \code{logratio} is extracted from \code{object}).
#'
#' Different prediction distances are proposed for discriminant analysis. The
#' reason is that our supervised models work with a dummy indicator matrix of
#' \code{Y} to indicate the class membership of each sample. The prediction of
#' a new observation results in either a predicted dummy variable (output
#' \code{object$predict}), or a predicted variate (output
#' \code{object$variates}). Therefore, an appropriate distance needs to be
#' applied to those predicted values to assign the predicted class. We propose
#' distances such as `maximum distance' for the predicted dummy variables,
#' `Mahalanobis distance' and `Centroids distance' for the predicted variates.
#'
#' \code{"max.dist"} is the simplest method to predict the class of a test
#' sample. For each new individual, the class with the largest predicted dummy
#' variable is the predicted class. This distance performs well in single data
#' set analysis with multiclass problems (PLS-DA).
#'
#' \code{"centroids.dist"} allocates to the new observation the class that
#' mimimises the distance between the predicted score and the centroids of the
#' classes calculated on the latent components or variates of the trained
#' model.
#'
#' \code{"mahalanobis.dist"} allocates the new sample the class defined as the
#' centroid distance, but using the Mahalanobis metric in the calculation of
#' the distance.
#'
#' In practice we found that the centroid-based distances
#' (\code{"centroids.dist"} and \code{"mahalanobis.dist"}), and specifically
#' the Mahalanobis distance led to more accurate predictions than the maximum
#' distance for complex classification problems and N-integration problems
#' (block.splsda). The centroid distances consider the prediction in
#' dimensional space spanned by the predicted variates, while the maximum
#' distance considers a single point estimate using the predicted scores on the
#' last dimension of the model. The user can assess the different distances,
#' and choose the prediction distance that leads to the best performance of the
#' model, as highlighted from the tune and perf outputs
#'
#' More (mathematical) details about the prediction distances are available in
#' the supplemental of the mixOmics article (Rohart et al 2017).
#'
#' For a visualisation of those prediction distances, see
#' \code{background.predict} that overlays the prediction area in
#' \code{plotIndiv} for a sPLS-DA object.
#'
#' %allocates the individual \eqn{x} to the class of \eqn{Y} minimizing
#' \eqn{dist(\code{x-variate}, G_l)}, where \eqn{G_l}, \eqn{l = 1,...,L} are
#' the centroids of the classes calculated on the \eqn{X}-variates of the
#' model. \code{"mahalanobis.dist"} allocates the individual \eqn{x} to the
#' class of \eqn{Y} as in \code{"centroids.dist"} but by using the Mahalanobis
#' metric in the calculation of the distance.
#'
#' For MINT objects, the \code{study.test} argument is required and provides
#' the grouping factor of \code{newdata}.
#'
#' For multi block analysis (thus block objects), \code{newdata} is a list of
#' matrices whose names are a subset of \code{names(object$X)} and missing
#' blocks are allowed. Several predictions are returned, either for each block
#' or for all blocks. For non discriminant analysis, the predicted values
#' (\code{predict}) are returned for each block and these values are combined
#' by average (\code{AveragedPredict}) or weighted average
#' (\code{WeightedPredict}), using the weights of the blocks that are
#' calculated as the correlation between a block's components and the outcome's
#' components.
#'
#' For discriminant analysis, the predicted class is returned for each block
#' (\code{class}) and each distance (\code{dist}) and these predictions are
#' combined by majority vote (\code{MajorityVote}) or weighted majority vote
#' (\code{WeightedVote}), using the weights of the blocks that are calculated
#' as the correlation between a block's components and the outcome's
#' components. NA means that there is no consensus among the block. For PLS-DA
#' and sPLS-DA objects, the prediction area can be visualised in plotIndiv via
#' the \code{background.predict} function.
#'
#' @aliases predict.pls predict.spls predict.plsda predict.splsda
#' predict.mint.pls predict.mint.spls predict.mint.plsda predict.mint.splsda
#' predict.mint.block.pls predict.mint.block.spls predict.mint.block.plsda
#' predict.mint.block.splsda
#' @param object object of class inheriting from
#' \code{"(mint).(block).(s)pls(da)"}.
#' @param newdata data matrix in which to look for for explanatory variables to
#' be used for prediction. Please note that this method does not perform
#' multilevel decomposition or log ratio transformations, which need to be
#' processed beforehand.
#' @param study.test For MINT objects, grouping factor indicating which samples
#' of \code{newdata} are from the same study. Overlap with \code{object$study}
#' are allowed.
#' @param dist distance to be applied for discriminant methods to predict the
#' class of new data, should be a subset of \code{"centroids.dist"},
#' \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details). Defaults to
#' \code{"all"}.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements). A numeric matrix or data frame. For a one level factor
#' decomposition, the input is a vector indicating the repeated measures on
#' each individual, i.e. the individuals ID. For a two level decomposition with
#' splsda models, the two factors are included in Y. Finally for a two level
#' decomposition with spls models, 2nd AND 3rd columns in design indicate those
#' factors (see example in \code{?splsda} and \code{?spls}).
#' @param ... not used currently.
#' @return \code{predict} produces a list with the following components:
#' \item{predict}{predicted response values. The dimensions correspond to the
#' observations, the response variables and the model dimension, respectively.
#' For a supervised model, it corresponds to the predicted dummy variables.}
#' \item{variates}{matrix of predicted variates.} \item{B.hat}{matrix of
#' regression coefficients (without the intercept).}
#'
#' \item{AveragedPredict}{if more than one block, returns the average predicted
#' values over the blocks (using the \code{predict} output)}
#' \item{WeightedPredict}{if more than one block, returns the weighted average
#' of the predicted values over the blocks (using the \code{predict} and
#' \code{weights} outputs)}
#'
#' \item{class}{predicted class of \code{newdata} for each
#' \eqn{1,...,}\code{ncomp} components.}
#'
#' \item{MajorityVote}{if more than one block, returns the majority class over
#' the blocks. NA for a sample means that there is no consensus on the
#' predicted class for this particular sample over the blocks.}
#' \item{WeightedVote}{if more than one block, returns the weighted majority
#' class over the blocks. NA for a sample means that there is no consensus on
#' the predicted class for this particular sample over the blocks.}
#' \item{weights}{Returns the weights of each block used for the weighted
#' predictions, for each nrepeat and each fold}
#'
#' \item{centroids}{matrix of coordinates for centroids.} \item{dist}{type of
#' distance requested.} \item{vote}{majority vote result for multi block
#' analysis (see details above).}
#' @author Florian Rohart, Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{mint.pls}}, \code{\link{mint.spls}},
#' \code{\link{mint.plsda}}, \code{\link{mint.splsda}},
#' \code{\link{block.pls}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}}, \code{\link{block.splsda}},
#' \code{\link{mint.block.pls}}, \code{\link{mint.block.spls}},
#' \code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}} and
#' visualisation with \code{\link{background.predict}} and
#' http://www.mixOmics.org for more details.
#' @references
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#'
#' Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris:
#' Editions Technic.
#' @keywords regression multivariate
#' @examples
#'
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")
#'
#' indiv1 <- c(200, 40, 60)
#' indiv2 <- c(190, 45, 45)
#' newdata <- rbind(indiv1, indiv2)
#' colnames(newdata) <- colnames(X)
#' newdata
#'
#' pred <- predict(linn.pls, newdata)
#'
#' plotIndiv(linn.pls, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE)
#' points(pred$variates[, 1], pred$variates[, 2], pch = 19, cex = 1.2)
#' text(pred$variates[, 1], pred$variates[, 2],
#' c("new ind.1", "new ind.2"), pos = 3)
#'
#' ## First example with plsda
#' X <- liver.toxicity$gene
#' Y <- as.factor(liver.toxicity$treatment[, 4])
#'
#'
#' ## if training is perfomed on 4/5th of the original data
#' samp <- sample(1:5, nrow(X), replace = TRUE)
#' test <- which(samp == 1)   # testing on the first fold
#' train <- setdiff(1:nrow(X), test)
#'
#' plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
#' test.predict <- predict(plsda.train, X[test, ], dist = "max.dist")
#' Prediction <- test.predict$class$max.dist[, 2]
#' cbind(Y = as.character(Y[test]), Prediction)
#'
#' \dontrun{
#' ## Second example with splsda
#' splsda.train <- splsda(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
#' test.predict <- predict(splsda.train, X[test, ], dist = "max.dist")
#' Prediction <- test.predict$class$max.dist[, 2]
#' cbind(Y = as.character(Y[test]), Prediction)
#'
#'
#' ## example with block.splsda=diablo=sgccda and a missing block
#' # need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
#' Y.mat = unmap(nutrimouse$diet)
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y.mat)
#' # with this design, all blocks are connected
#' design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
#' byrow = TRUE, dimnames = list(names(data), names(data)))
#'
#' # train on 75% of the data
#' ind.train=NULL
#' for(i in 1:nlevels(nutrimouse$diet))
#' ind.train=c(ind.train,which(nutrimouse$diet==levels(nutrimouse$diet)[i])[1:6])
#'
#' #training set
#' gene.train=nutrimouse$gene[ind.train,]
#' lipid.train=nutrimouse$lipid[ind.train,]
#' Y.mat.train=Y.mat[ind.train,]
#' Y.train=nutrimouse$diet[ind.train]
#' data.train=list(gene=gene.train,lipid=lipid.train,Y=Y.mat.train)
#'
#' #test set
#' gene.test=nutrimouse$gene[-ind.train,]
#' lipid.test=nutrimouse$lipid[-ind.train,]
#' Y.mat.test=Y.mat[-ind.train,]
#' Y.test=nutrimouse$diet[-ind.train]
#' data.test=list(gene=gene.test,lipid=lipid.test)
#'
#' # example with block.splsda=diablo=sgccda and a missing block
#' res.train = block.splsda(X=list(gene=gene.train,lipid=lipid.train),Y=Y.train,
#' ncomp=3,keepX=list(gene=c(10,10,10),lipid=c(5,5,5)))
#' test.predict = predict(res.train, newdata=data.test[2], method = "max.dist")
#'
#'
#' ## example with mint.splsda
#'
#' #training set
#' ind.test = which(stemcells$study == "3")
#' gene.train = stemcells$gene[-ind.test,]
#' Y.train = stemcells$celltype[-ind.test]
#' study.train = factor(stemcells$study[-ind.test])
#'
#' #test set
#' gene.test = stemcells$gene[ind.test,]
#' Y.test = stemcells$celltype[ind.test]
#' study.test = factor(stemcells$study[ind.test])
#'
#' res = mint.splsda(X = gene.train, Y = Y.train, ncomp = 3, keepX = c(10, 5, 15),
#' study = study.train)
#'
#' pred = predict(res, newdata = gene.test, study.test = study.test)
#'
#' data.frame(Truth = Y.test, prediction = pred$class$max.dist)
#'
#' }

#' @importFrom stats predict
#'
perf.mint.spls  = perf.mint.pls = function(object,
validation = c("Mfold", "loo"),
folds = 10,
progressBar = TRUE,
...)
{
    stop("Yet to be implemented")
}

# ---------------------------------------------------
# perf for mint.plsda and mint.splsda object
# ---------------------------------------------------
perf.mint.splsda = perf.mint.plsda = function (object,
dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
auc = FALSE,
progressBar = TRUE,
...
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#

    #------------------#
    #-- check entries --#
    X = object$X
    Y = object$Y
    study = object$study
    ncomp = object$ncomp
    scale = object$scale

    keepX = apply(object$loadings$X, 2, function(x){sum(x!=0)})

    tol = object$tol
    max.iter = object$max.iter

    dist = match.arg(dist, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    if (any(dist == "all"))
    {
        nmthdd = 3
        dist = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
        nmthdd = length(dist)
    }

    if (!is.logical(progressBar))
    stop("'progressBar' must be either TRUE or FALSE")

    #-- end checking --#
    #------------------#

    near.zero.var = !is.null(object$nzv) # if near.zero.var was used, we set it to TRUE. if not used, object$nzv is NULL


    # -------------------------------------
    # added: first check for near zero var on the whole data set
    if (near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]

            if (ncol(X)==0)
            stop("No more predictors after Near Zero Var has been applied!")
        }
    }
    # and then we start from the X data set with the nzv removed

    prediction.all = class.all = list()
    for(ijk in dist)
    {
        class.all[[ijk]] = matrix(nrow = nrow(X), ncol = ncomp,
        dimnames = list(rownames(X), c(paste0('comp', 1 : ncomp))))
    }

    if(auc)
    {
        auc.mean=list()
        auc.mean.study=list()
    }

    study.specific = global = list()
    for (study_i in 1:nlevels(study)) #LOO on the study factor
    {
        study.specific[[study_i]] =list()
        study.specific[[study_i]]$BER = global$BER = matrix(0,nrow = ncomp, ncol = length(dist),
        dimnames = list(c(paste0('comp', 1 : ncomp)), dist))

        study.specific[[study_i]]$overall = global$overall = matrix(0,nrow = ncomp, ncol = length(dist),
        dimnames = list(c(paste0('comp', 1 : ncomp)), dist))

        study.specific[[study_i]]$error.rate.class = list()
        for(ijk in dist)
        study.specific[[study_i]]$error.rate.class[[ijk]] = global$error.rate.class[[ijk]] = matrix(0,nrow = nlevels(Y), ncol = ncomp,
        dimnames = list(levels(Y),c(paste0('comp', 1 : ncomp))))

    }
    names(study.specific) =levels(study)

    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1 : ncomp)
    {

        already.tested.X = keepX[1:comp]


        if (progressBar == TRUE)
        cat("\ncomp",comp, "\n")

        #-- set up a progress bar --#
        if (progressBar ==  TRUE)
        {
            pb = txtProgressBar(style = 3)
            nBar = 1
        } else {
            pb = FALSE
        }

        M = nlevels(study)
        names.study = levels(study)
        features = NULL
        prediction.comp = matrix(0, nrow = nrow(X), ncol = nlevels(Y), dimnames = list(rownames(X), levels(Y)))

        class.comp = list()
        for(ijk in dist)
        class.comp[[ijk]] = matrix(0, nrow = nrow(X), ncol = 1)# prediction of all samples for each test.keepX and  nrep at comp fixed

        if(auc)
        auc.mean.study[[comp]] = list()

        for (study_i in 1:M) #LOO on the study factor
        {
            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (study_i-1)/M)

            omit = which(study %in% names.study[study_i])
            X.train = X[-omit,]
            Y.train = factor(Y[-omit])
            study.learn.CV = factor(as.character(study[-omit]))

            X.test = X[omit, , drop = FALSE] #note: drop is useless as there should always be more than a single sample in a study
            Y.test = Y[omit]
            study.test.CV = factor(as.character(study[omit]))

            #---------------------------------------#
            #-- near.zero.var ----------------------#
            if(near.zero.var == TRUE)
            {
                remove.zero = nearZeroVar(X.train)$Position

                if (length(remove.zero) > 0)
                {
                    X.train = X.train[, -c(remove.zero),drop = FALSE]
                    X.test = X.test[, -c(remove.zero),drop = FALSE]
                }
            }
            #-- near.zero.var ----------------------#
            #---------------------------------------#

            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (study_i-1)/M)

            object.res = mint.splsda(X.train, Y.train, study = study.learn.CV, ncomp = comp,
            keepX = already.tested.X,
            scale = scale, mode = "regression")

            test.predict.sw <- predict.mixo_spls(object.res, newdata = X.test, dist = dist, study.test = study.test.CV)
            prediction.comp[omit, match(levels(Y.train),levels(Y))] =  test.predict.sw$predict[, , comp]

            for(ijk in dist)
            class.comp[[ijk]][omit,1] =  test.predict.sw$class[[ijk]][, comp] #levels(Y)[test.predict.sw$class[[ijk]][, ncomp]]


            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (study_i)/M)

            # result per study
            #BER
            study.specific[[study_i]]$BER[comp,] = sapply(test.predict.sw$class, function(x){
            conf = get.confusion_matrix(truth = Y[omit], all.levels = levels(Y), predicted = x[,comp])
            get.BER(conf)
            })

            #overall
            study.specific[[study_i]]$overall[comp,] = sapply(test.predict.sw$class, function(x){
                conf = get.confusion_matrix(truth = Y[omit], all.levels = levels(Y), predicted = x[,comp])
                out = sum(apply(conf, 1, sum) - diag(conf)) / length(Y[omit])
            })

            #classification for each level of Y
            temp = lapply(test.predict.sw$class, function(x){
                conf = get.confusion_matrix(truth = Y[omit], all.levels = levels(Y), predicted = x[,comp])
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y[omit])
            })
            for (ijk in dist)
            study.specific[[study_i]]$error.rate.class[[ijk]][,comp] = temp[[ijk]]

            #AUC per study
            if(auc)
            {
                data = list()
                data$outcome = Y[omit]
                data$data = prediction.comp[omit, ]
                auc.mean.study[[comp]][[study_i]] = statauc(data)
            }

        } # end study_i 1:M (M folds)

        for (ijk in dist)
        {
            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[ijk]][,comp] = class.comp[[ijk]][,1]
        }
        prediction.all[[comp]] = prediction.comp


        # global results
        #BER
        global$BER[comp,] = sapply(class.comp, function(x){
            conf = get.confusion_matrix(truth = factor(Y), predicted = x)
            get.BER(conf)
        })

        #overall
        global$overall[comp,] = sapply(class.comp, function(x){
            conf = get.confusion_matrix(truth = factor(Y), predicted = x)
            out = sum(apply(conf, 1, sum) - diag(conf)) / length(Y)
        })

        #classification for each level of Y
        temp = lapply(class.comp, function(x){
            conf = get.confusion_matrix(truth = factor(Y), predicted = x)
            out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
        })
        for (ijk in dist)
        global$error.rate.class[[ijk]][,comp] = temp[[ijk]]

        #AUC global
        if(auc)
        {
            names(auc.mean.study[[comp]]) = names.study

            data = list()
            data$outcome = Y
            data$data = prediction.comp
            auc.mean[[comp]] = statauc(data)
        }


    } # end comp
    names(prediction.all) = paste0('comp', 1:ncomp)

    result = list(study.specific.error = study.specific,
    global.error = global,


    predict = prediction.all,
    class = class.all)

    if(auc)
    {
        names(auc.mean) = names(auc.mean.study) = paste0('comp', 1:ncomp)
        result$auc = auc.mean
        result$auc.study = auc.mean.study
    }

    if (progressBar == TRUE)
    cat('\n')


    # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
    if(nlevels(study) > 2 & ncomp >1)
    {
        measure = c("overall","BER")
        ncomp_opt = matrix(NA, nrow = length(measure), ncol = length(dist),
        dimnames = list(measure, dist))

        for (measure_i in measure)
        {
            for (ijk in dist)
            {
                mat.error.rate = sapply(study.specific, function(x){x[[measure_i]][,ijk]})
                ncomp_opt[measure_i, ijk] = t.test.process(t(mat.error.rate))
            }
        }
    } else {
        ncomp_opt = NULL
    }

    result$choice.ncomp = ncomp_opt


    # added
    if (near.zero.var == TRUE)
    result$nzvX = nzv$Position

    class(result) = c("perf","perf.mint.splsda.mthd")
    result$call = match.call()

    return(invisible(result))
}
