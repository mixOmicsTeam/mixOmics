#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

#
# created: 2013
# last modified: 22-04-2016
#
# Copyright (C) 2013
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







#' Tuning functions for multilevel sPLS method
#'
#' For a multilevel spls analysis, the tuning criterion is based on the
#' maximisation of the correlation between the components from both data sets
#'
#' For a multilevel spls analysis, the tuning criterion is based on the
#' maximisation of the correlation between the components from both data sets
#'
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y \code{if(method = 'spls')} numeric vector or matrix of continuous
#' responses (for multi-response models) \code{NA}s are allowed.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements) that indicates the repeated measures on each individual, i.e.
#' the individuals ID. See Details.
#' @param ncomp the number of components to include in the model.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}.
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set
#' @param test.keepY numeric vector for the different number of variables to
#' test from the \eqn{Y} data set
#' @param already.tested.X Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{X} data set on
#' the firsts components.
#' @param already.tested.Y Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{Y} data set on
#' the firsts components.
#' @return
#'
#' \item{cor.value}{correlation between latent variables}
#' @author Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart.
#' @seealso \code{\link{splsda}}, \code{\link{predict.splsda}} and
#' http://www.mixOmics.org for more details.
#' @references mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @examples
#'
#' # note: we made up those data, pretending they are repeated measurements
#' repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
#' 6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
#' 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
#' 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
#' summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
#'
#' # this is a spls (unsupervised analysis) so no need to mention any factor in design
#' # we only perform a one level variation split
#' design <- data.frame(sample = repeat.indiv)
#'
#' tune.splslevel(X = liver.toxicity$gene,
#' Y=liver.toxicity$clinic,
#' multilevel = design,
#' test.keepX = c(5,10,15),
#' test.keepY = c(1,2,5),
#' ncomp = 1)
#'
#'
#' @export tune.splslevel
tune.splslevel <- function (X, Y,
multilevel,
ncomp = NULL,
mode="regression",
test.keepX = rep(ncol(X), ncomp),
test.keepY = rep(ncol(Y), ncomp),
already.tested.X = NULL,
already.tested.Y = NULL)
{

    message("For a multilevel spls analysis, the tuning criterion is based on the maximisation of the correlation between the components from both data sets")

    Y = as.matrix(Y)
    if (length(dim(Y)) != 2 || !is.numeric(Y))
    stop("'Y' must be a numeric matrix.")

    if (!is.null(already.tested.X))
    cat("Number of X variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X, "\n")

    if (!is.null(already.tested.Y))
    cat("Number of Y variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.Y, "\n")

    if ((!is.null(already.tested.X)) && is.null(already.tested.Y))
    stop("Input already.tested.Y is missing")

    if ((!is.null(already.tested.Y)) && is.null(already.tested.X))
    stop("Input already.tested.X is missing")

    if (length(already.tested.X) != (ncomp - 1))
    stop("The number of already.tested.X parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)

    if (length(already.tested.Y) != (ncomp - 1))
    stop("The number of already.tested.Y parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)

    if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X)))
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)

    if ((!is.null(already.tested.Y)) && (!is.numeric(already.tested.Y)))
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)

    Xw <- suppressMessages(withinVariation(X = X, design = multilevel))
    Yw <- suppressMessages(withinVariation(X = Y, design = multilevel))

    cor.value = matrix(nrow = length(test.keepX), ncol = length(test.keepY))
    rownames(cor.value) = paste("varX ", test.keepX, sep = "")
    colnames(cor.value) = paste("varY ", test.keepY, sep = "")

    for (i in 1:length(test.keepX))
    {
        for (j in 1:length(test.keepY))
        {
            if (ncomp == 1)
            {
                spls.train = mixOmics::spls(Xw, Yw, ncomp = ncomp,
                keepX = test.keepX[i],
                keepY = test.keepY[j],
                mode = mode)
            } else {
                spls.train = mixOmics::spls(Xw, Yw, ncomp = ncomp,
                keepX = c(already.tested.X, test.keepX[i]),
                keepY = c(already.tested.Y, test.keepY[j]),
                mode = mode)
            }

            cor.value[i, j] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
        }
    }
    return(list(cor.value = cor.value))
}
