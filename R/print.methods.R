#############################################################################################################
# Author :
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Leigh Coonan, Student, University of Queensland, Australia
#   Fangzhou Yao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia and
#   Shangai University of Finance and Economics, Shanghai, P.R. China
#   Jeff Coquery, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia and
#   Sup Biotech, Paris, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2009
# last modified: 05-10-2017
#
# Copyright (C) 2009
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
#' Print Methods for CCA, (s)PLS, PCA and Summary objects
#'
#' Produce \code{print} methods for class \code{"rcc"}, \code{"pls"},
#' \code{"spls"}, \code{"pca"}, \code{"rgcca"}, \code{"sgcca"} and
#' \code{"summary"}.
#'
#' \code{print} method for \code{"rcc"}, \code{"pls"}, \code{"spls"}
#' \code{"pca"}, \code{"rgcca"}, \code{"sgcca"} class, returns a description of
#' the \code{x} object including: the function used, the regularization
#' parameters (if \code{x} of class \code{"rcc"}), the (s)PLS algorithm used
#' (if \code{x} of class \code{"pls"} or \code{"spls"}), the samples size, the
#' number of variables selected on each of the sPLS components (if \code{x} of
#' class \code{"spls"}) and the available components of the object.
#'
#' \code{print} method for \code{"summary"} class, gives the (s)PLS algorithm
#' used (if \code{x} of class \code{"pls"} or \code{"spls"}), the number of
#' variates considered, the canonical correlations (if \code{x} of class
#' \code{"rcc"}), the number of variables selected on each of the sPLS
#' components (if \code{x} of class \code{"spls"}) and the available components
#' for Communalities Analysis, Redundancy Analysis and Variable Importance in
#' the Projection (VIP).
#'
#' @aliases print print.rcc print.mixo_pls print.mixo_spls print.summary
#' print.pca print.spca print.rgcca print.sgcca
#' @param x object of class inherited from \code{"rcc"}, \code{"pls"},
#' \code{"spls"}, \code{"pca"}, \code{"spca"}, \code{"rgcca"}, \code{"sgcca"}or
#' \code{"summary"}.
#' @param list() not used currently.
#' @return none
#' @author Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.
#' @seealso \code{\link{rcc}}, \code{\link{pls}}, \code{\link{spls}},
#' \code{\link{vip}}.
#' @keywords regression multivariate
#' @examples
#'
#' ## print for objects of class 'rcc'
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
#' print(nutri.res)
#'
#' \dontrun{
#' ## print for objects of class 'summary'
#' more <- summary(nutri.res, cutoff = 0.65)
#' print(more)
#'
#' ## print for objects of class 'pls'
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y)
#' print(linn.pls)
#'
#' ## print for objects of class 'spls'
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
#' keepY = c(10, 10, 10))
#' print(toxicity.spls)
#' }
#'

#------------------ print method for pls ------------------#
print.mixo_pls <-
function(x, ...)
{

    mode = paste("'", x$mode, "'", sep = "")



    cat(" PLS with a", mode, "mode with", x$ncomp, "PLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")

    cat(" No variable selection. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")

}

#------------------ print method for mint.pls ------------------#
print.mint.pls <-
function(x, ...)
{

    mode = paste("'", x$mode, "'", sep = "")



    cat(" MINT PLS with a", mode, "mode with", x$ncomp, "PLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
    cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")

    cat(" No variable selection. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" loading vectors per study: see object$loadings.partial \n")
    cat(" variates: see object$variates \n")
    cat(" variates per study: see object$variates.partial \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
}

#------------------ print method for plsda ------------------#
print.mixo_plsda <-
function(x, ...)
{



    cat(" PLS-DA (regression mode) with", x$ncomp, "PLS-DA components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")

    cat(" No variable selection. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow, cim \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" auc \n")

}

#------------------ print method for mint.plsda ------------------#
print.mint.plsda <-
function(x, ...)
{



    cat(" MINT PLS-DA (regression mode) with", x$ncomp, "PLS-DA components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
    cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")

    cat(" No variable selection. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" loading vectors per study: see object$loadings.partial \n")
    cat(" variates: see object$variates \n")
    cat(" variates per study: see object$variates.partial \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow, cim \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" perf, auc\n")
}

#----------------- print method for spls ------------------#
print.mixo_spls <-
function(x, ...)
{

    mode = paste("'", x$mode, "'", sep = "")
    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")



    cat(" sPLS with a", mode, "mode with", x$ncomp, "sPLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")

    cat(" Selection of", keepX, "variables on each of the sPLS components on the X data set. \n")
    cat(" Selection of", keepY, "variables on each of the sPLS components on the Y data set. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
}

#----------------- print method for mint.spls ------------------#
print.mint.spls <-
function(x, ...)
{

    mode = paste("'", x$mode, "'", sep = "")
    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")



    cat(" MINT sPLS with a", mode, "mode with", x$ncomp, "sPLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
    cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")

    cat(" Selection of", keepX, "variables on each of the sPLS components on the X data set. \n")
    cat(" Selection of", keepY, "variables on each of the sPLS components on the Y data set. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" loading vectors per study: see object$loadings.partial \n")
    cat(" variates: see object$variates \n")
    cat(" variates per study: see object$variates.partial \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar\n")
}


#----------------- print method for splsda ------------------#
print.mixo_splsda <-
function(x, ...)
{

    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")


    cat(" sPLS-DA (regression mode) with", x$ncomp, "sPLS-DA components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")

    cat(" Selection of", keepX, "variables on each of the sPLS-DA components on the X data set. \n")
    cat(" No Y variables can be selected. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow, cim \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar, tune, perf, auc \n")
}

#----------------- print method for mint.splsda ------------------#
print.mint.splsda <-
function(x, ...)
{

    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")


    cat(" MINT sPLS-DA (regression mode) with", x$ncomp, "sPLS-DA components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
    cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")

    cat(" Selection of", keepX, "variables on each of the sPLS-DA components on the X data set. \n")
    cat(" No Y variables can be selected. \n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" loading vectors per study: see object$loadings.partial \n")
    cat(" variates: see object$variates \n")
    cat(" variates per study: see object$variates.partial \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow, cim \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, cim \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar, tune, perf, auc \n")

}

#------------------ print method for rcc ------------------#
print.rcc <-
function(x, ...)
{



    cat(" rCCA with", x$ncomp, "components and regularization parameters", x$lambda[1], "and", x$lambda[2], "for the X and Y data. \n")
    cat(" You entered data X of dimensions :", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions :", nrow(x$Y), ncol(x$Y), "\n\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" canonical correlations: see object$cor \n")
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}

# ------------------------ print for pca --------------------------------
print.pca <- function(x, ...)
{

    ind.show = min(10, x$ncomp)

    cat("Eigenvalues for the first", ind.show, "principal components, see object$sdev^2:", "\n")
    print((x$sdev[1:ind.show])^2)
    cat("\n")

    per.var = x$explained_variance
    cum.var=as.vector(cumsum(per.var))
    x$sdev=as.vector(x$sdev)
    names(x$sdev) = paste("PC", 1:length(x$sdev), sep = "")
    names(per.var) = paste("PC", 1:length(per.var), sep = "")
    names(cum.var) = paste("PC", 1:length(cum.var), sep = "")

    cat("Proportion of explained variance for the first", ind.show, "principal components, see object$explained_variance:", "\n")
    print(per.var[1:ind.show])
    cat("\n")

    cat("Cumulative proportion explained variance for the first", ind.show, "principal components, see object$cum.var:", "\n")
    print(cum.var[1:ind.show])
    cat("\n")

    cat(" Other available components: \n", "-------------------- \n")
    cat(" loading vectors: see object$rotation \n")
}

# ------------------------ print for spca -------------------------
print.spca <-
function(x, ...)
{
    cat(" sparse PCA with", x$ncomp, "principal components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")

    cat(" Selection of", x$keepX, "variables on each of the principal components on the X data set. \n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" loading vectors: see object$rotation \n")
    cat(" principal components: see object$x \n")
    cat(" cumulative explained variance: see object$varX \n")
    cat(" variable names: see object$names \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar, tune\n")

}

# ------------------------ print for ipca -------------------------
print.ipca <-
function(x, ...)
{
    cat(" IPCA with", x$ncomp, "independent components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" unmixing matrix: see object$unmixing \n")
    cat(" independent principal components: see object$x \n")
    cat(" mxing matrix: see object$mixing \n")
    cat(" kurtosis: see object$kurtosis \n")
    cat(" variable names: see object$names \n")
    cat(" independent loading vectors: see object$loadings \n")
}

# ------------------------ print for sipca -------------------------
print.sipca <-
function(x, ...)
{

    cat(" Sparse IPCA with", x$ncomp, "independent components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")

    cat(" Selection of", x$keepX, "variables on each of the principal components on the X data set. \n")

    cat(" Main numerical outputs: \n",
    "-------------------- \n")

    cat(" unmixing matrix: see object$unmixing \n")
    cat(" independent principal components: see object$x \n")
    cat(" mxing matrix: see object$mixing \n")
    cat(" kurtosis: see object$kurtosis \n")
    cat(" variable names: see object$names \n")
    cat(" independent loading vectors: see object$loadings \n")

}

# ------------------------ print for rgcca -------------------------
print.rgcca <- function(x, ...)
{

    # components
    for(k in 1:length(x$blocks)){
        cat(" rGCCA with", x$ncomp[k], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat("\n")

    # dimension
    for(k in 1 : length(x$blocks)){
        cat(" Dimension of block", k, 'is ', dim(x$blocks[[k]]), "\n")
    }
    cat("\n")
    cat(" Main numerical outputs: \n", "-------------------- \n")

    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")

    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar\n")

}


# ------------------------ print for sgcca -------------------------
print.sgcca<- function(x, ...)
{

    # components
    for(k in 1 : length(x$X)){
        cat(" sGCCA with", x$ncomp[[k]], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat("\n")

    # dimension
    for(k in 1 : length(x$X)){
        cat(" Dimension of block", k, 'is ', dim(x$X[[k]]), "\n")
    }
    cat("\n")

    # selected variables
    list.select = list()
    for(k in 1:length(x$X)){
        list.select[[k]] = apply(x$loadings[[k]], 2, function(x){sum(x!=0)})
        cat(" Selection of", list.select[[k]], "variables on each of the sGCCA components on the block", k, "\n")
    }
    cat("\n")
    cat(" Main numerical outputs: \n", "-------------------- \n")
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")


    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network\n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar \n")

}


# ------------------------ print for sgcca -------------------------
print.sgccda<- function(x, ...)
{

    # components
    for(k in 1 : length(x$X)){
        cat(" sGCCA with", x$ncomp[[k]], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat(" sGCCA with", x$ncomp[x$indY], "components on the outcome Y\n")
    cat("\n")


    # dimension
    for(k in 1 : length(x$X)){
        cat(" Dimension of block", k, 'is ', dim(x$X[[k]]), "\n")
    }
    cat(" Outcome Y has", nlevels(x$Y), "levels \n")
    cat("\n")

    # selected variables
    list.select = list()
    for(k in 1:length(x$X)){
        list.select[[k]] = apply(x$loadings[[k]], 2, function(x){sum(x!=0)})
        cat(" Selection of", list.select[[k]], "variables on each of the sGCCA components on the block", k, "\n")
    }
    cat("\n")
    cat(" Main numerical outputs: \n", "-------------------- \n")
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")


    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow, cimDiablo, plotDiablo \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, circosPlot \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar, perf, auc \n")

}


#------- print for summary with (s)PLS object or rcc ---------#
print.summary <-
function(x, ...)
{

    print.gap = 4
    what = x$what
    digits = x$digits

    #--------------------- output pls/spls ---------------------#
    if(is(x, c("pls", "spls"))){

        if (is(x, "pls"))
        {
            cat(" PLS mode:", x$mode)
            cat("\n Number of variates considered:", x$ncomp, "\n")
        } else {
            cat(" sPLS mode:", x$mode)
            cat("\n Number of variates considered:", x$ncomp)
            cat("\n Number of X-variables selected on each of the sPLS components:", x$keepX)
            cat("\n Number of Y-variables selected on each of the sPLS components:", x$keepY, "\n")
        }

        #---------- affichage communaute ----------#
        if (any(what == "all") || any(what == "communalities"))
        {
            cat("\n\n Communalities Analysis:\n",
            "----------------------")

            cat("\n X-Variables vs their own Variates: see object$CM.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$CM.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$CM.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$CM.Y$opp \n")
        }

        #--------- affichage redondance -----------#
        if (any(what == "all") || any(what == "redundancy"))
        {
            cat("\n\n Redundancy Analysis:\n",
            "-------------------\n")

            cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
        }

        #---------- tableau VIP ---------#
        if (any(what == "all") || any(what == "VIP"))
        {
            cat("\n\n", "Variable Importance in the Projection (VIP): see object$VIP \n",
            "------------------------------------------- \n\n")
        }

    }  #end if pls

    # ---------------------- output rcc ------------------------#
    if(is(x, "rcc"))
    {
        print.gap = 4
        if (any(what == "all"))
        {
            cat(" Number of canonical variates considered:", x$ncomp, "\n")
            cat("\n Canonical correlations:",
            "\n ----------------------\n")
            print(round(x$can.cor, digits = digits), print.gap = print.gap)
        }

        #-- affichage communaute --#
        if (any(what == "all") || any(what == "communalities"))
        {
            cat("\n\n Canonical Communalities Analysis:\n",
            "--------------------------------")

            cat("\n X-Variables vs their own Canonical Variates: see object$Cm.X$own \n")
            cat("\n X-Variables vs the opposite Canonical Variates: see object$Cm.X$opp \n")
            cat("\n Y-Variables vs their own Canonical Variates: see object$Cm.Y$own \n")
            cat("\n Y-Variables vs the opposite Canonical Variates: see object$Cm.Y$opp \n")
        }

        #--------- affichage redondance -----------#
        if (any(what == "all") || any(what == "redundancy"))
        {
            cat("\n\n Redundancy Analysis:\n",
            "-------------------\n")

            cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
        }

    }  #end rcc
}


# perf.diablo / sgccda.mthd
# perf.splsda = perf.plsda / splsda.mthd plsda.mthd
# perf.spls  = perf.pls / spls.mthd pls.mthd

print.perf.pls.mthd = function(x, ...)
{
    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" MSEP, R2, Q2, Q2.total, RSS, PRESS. See the help file ?perf \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")

}

print.perf.spls.mthd = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" MSEP, R2, Q2, Q2.total, RSS, PRESS. See the help file ?perf \n")
    cat(" Stable features of X on each component: see object$features$stable.X \n")
    cat(" Stable features of Y on each component: see object$features$stable.Y \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")

}


print.perf.plsda.mthd = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Error rate (overall or BER) for each component and for each distance: see object$error.rate \n")
    cat(" Error rate per class, for each component and for each distance: see object$error.rate.class \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" AUC values: see object$auc if auc = TRUE \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")

}

print.perf.splsda.mthd = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Error rate (overall or BER) for each component and for each distance (averaged over the nrepeat): see object$error.rate \n")
    cat(" Error rate (overall or BER) for each component, for each distance and for each repeat: see object$error.rate.all \n")
    cat(" Error rate per class, for each component and for each distance: see object$error.rate.class \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" Stable features on each component: see object$features$stable \n")
    cat(" AUC values: see object$auc if auc = TRUE \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")

}

print.perf.mint.splsda.mthd = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Study-specific error rate (overall, BER and error rate per class) for each component and for each distance: see object$study.specific.error \n")
    cat(" Global error rate (overall, BER and error rate per class) for each component and for each distance: see object$global.error \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" AUC values: see object$auc \n")
    cat(" AUC values per study: see object$auc.study if auc = TRUE \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")

}

print.perf.sgccda.mthd = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Error rate (overall or BER) for each component and for each distance: see object$error.rate \n")
    cat(" Error rate per class, for each component and for each distance: see object$error.rate.class \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" Stable features on each component: see object$features$stable \n")
    cat(" AUC values: see object$auc if auc = TRUE \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")

}


# tune: "spls", "splsda", "mint.splsda", "rcc", "pca"
print.tune.pca = function(x, ...)
{

    cat(" for all principal components, see object$sdev, object$explained_variance and object$cum.var\n")
}

print.tune.rcc = function(x, ...)
{

    cat("  lambda1 = ", x$opt.lambda1, ", see object$opt.lambda1\n", " lambda2 = ", x$opt.lambda2, ",  see object$opt.lambda2\n",
    "CV-score = ", x$opt.score, ", see object$opt.score\n")
}

print.tune.splsda = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n")
    cat(" Optimal number of components: see object$choice.ncomp \n")
    cat(" Error rate for each tested keepX and for each component (averaged over the nrepeat, mean and standard deviation): see object$error.rate and object$error.rate.sd  \n")
    cat(" Error rate for each tested keepX, for each component and for each repeat: see object$error.rate.all \n")
    cat(" Error rate per class obtained with the optimal keepX, for each component and for each nrepeat: see object$error.rate.class \n")
    cat(" AUC when applicable, (note that those results may differ with sPLS-DA prediction, see details): see object$AUC \n\n")

    cat(" Other outputs are available, and details on those outputs in ?tune.splsda.  \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
}

print.tune.mint.splsda = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n")
    cat(" Error rate for each tested keepX and for each component: see object$error.rate \n")
    cat(" Error rate per class obtained with the optimal keepX, for each component: see object$error.rate.class \n\n")

    cat(" Other outputs available, see ?tune.mint.splsda \n\n")

    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
}

print.tune.block.splsda = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n")
    cat(" Error rate for each tested keepX and for each component (averaged over the nrepeat): see object$error.rate \n")
    cat(" Error rate for each tested keepX, for each component and for each repeat: see object$error.rate.all \n")
    cat(" Error rate per class obtained with the optimal keepX, for each component and for each nrepeat: see object$error.rate.class \n\n")

    cat(" Other outputs available, see ?tune.splsda.  \n\n")

}

print.predict = function(x, ...)
{

    cat(" Main numerical outputs: \n",
    "-------------------- \n")
    if(is.list(x$predict)) #block analysis
    {
        cat(" Prediction values of the test samples for each block and each component: see object$predict \n")
        cat(" variates of the test samples for each block and each component: see object$variates \n")


        if(!is.null(x$dist)) #DA object, more outputs
        {
            cat(" Classification of the test samples for each distance, for each block and each component: see object$class \n")
            cat(" Majority vote of the test samples for each distance and each component: see object$vote \n")
        }
    }else{
        cat(" Prediction values of the test samples for each component: see object$predict \n")
        cat(" variates of the test samples: see object$variates \n")

        if(!is.null(x$dist)) #DA object, more outputs
        cat(" Classification of the test samples for each distance and for each component: see object$class \n")
    }

}


