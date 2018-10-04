#############################################################################################################
# Author :
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2009
# last modified:
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


#--------------------------------------------------------#
#-- Includes summary.pls, summary.spls and summary.rcc --#
#--------------------------------------------------------#

#--------------------- PLS and sPLS ---------------------#
summary.mixo_pls <- summary.mixo_spls <-
function(object, 
         what = c("all", "communalities", "redundancy", "VIP"), 
         digits = 4, 
         keep.var = FALSE,
	 ...) 
{

    #-- initialisation des matrices --#
    ncomp = object$ncomp
    n = nrow(object$X)
    p = ncol(object$X)
    q = ncol(object$Y)
    result = list( )
     
    if (keep.var) {
        gp.X = apply(object$loadings$X, 1, sum) != 0
        gp.Y = apply(object$loadings$Y, 1, sum) != 0
    }
    else {
        gp.X = 1:p
        gp.Y = 1:q
    }
     
    #-- communaut? --#
    if (any(what == "all") || any(what == "communalities")) {
        # communaut? inter-groupe #
        #-------------------------#
        Cm.XvsV = cor(object$X[, gp.X], object$variate$Y[, 1:ncomp], use = "pairwise")
        Cm.XvsV = t(apply(Cm.XvsV^2, 1, cumsum))
        colnames(Cm.XvsV) = paste("comp", 1:ncomp)
         
        Cm.YvsU = cor(object$Y[, gp.Y], object$variate$X[, 1:ncomp], use = "pairwise")
        Cm.YvsU = t(apply(Cm.YvsU^2, 1, cumsum))
        colnames(Cm.YvsU) = paste("comp", 1:ncomp)
         
        # communaut? intra-groupe #
        #-------------------------#
        Cm.XvsU = cor(object$X[, gp.X], object$variates$X[, 1:ncomp], use = "pairwise")
        Cm.XvsU = t(apply(Cm.XvsU^2, 1, cumsum))
        colnames(Cm.XvsU) = paste("comp", 1:ncomp)
         
        Cm.YvsV = cor(object$Y[, gp.Y], object$variates$Y[, 1:ncomp], use = "pairwise")
        Cm.YvsV = t(apply(Cm.YvsV^2, 1, cumsum))
        colnames(Cm.YvsV) = paste("comp", 1:ncomp)
         
        result$Cm.X = list(own = Cm.XvsU, opp = Cm.XvsV)
        result$Cm.Y = list(own = Cm.YvsV, opp = Cm.YvsU)
    }
     
    #-- redondance --#
    if (any(what == "all") || any(what == "redundancy")) {
        Rd.XvsU = cor(object$X[, gp.X], object$variates$X[, 1:ncomp], use = "pairwise")
        Rd.XvsU = apply(Rd.XvsU^2, 2, sum)/p
        Rd.XvsV = cor(object$X[, gp.X], object$variates$Y[, 1:ncomp], use = "pairwise")
        Rd.XvsV = apply(Rd.XvsV^2, 2, sum)/p
         
        Rd.YvsU = cor(object$Y[, gp.Y], object$variates$X[, 1:ncomp], use = "pairwise")
        Rd.YvsU = apply(Rd.YvsU^2, 2, sum)/q
        Rd.YvsV = cor(object$Y[, gp.Y], object$variates$Y[, 1:ncomp], use = "pairwise")
        Rd.YvsV = apply(Rd.YvsV^2, 2, sum)/q
         
        own = cbind(Rd.XvsU, cumsum(Rd.XvsU))
        colnames(own) = c("Proportion", "Cumulative")
        opp = cbind(Rd.XvsV, cumsum(Rd.XvsV))
        colnames(opp) = c("Proportion", "Cumulative")
        result$Rd.X = list(own = own, opp = opp)
         
        own = cbind(Rd.YvsV, cumsum(Rd.YvsV))
        colnames(own) = c("Proportion", "Cumulative")
        opp = cbind(Rd.YvsU, cumsum(Rd.YvsU))
        colnames(opp) = c("Proportion", "Cumulative")
        result$Rd.Y = list(own = own, opp = opp)
    }
     
    #-- affichage --#
    if (any(what == "all") || any(what == "communalities") || 
        any(what == "redundancy") || any(what == "VIP")) {
        result$ncomp = ncomp
        result$mode = object$mode
        result$keep.var = list(X = colnames(object$X)[gp.X], Y = colnames(object$Y)[gp.Y])
    }
     
    #-- tableau VIP --#
    if (any(what == "all") || any(what == "VIP")) {
        VIP = vip(object)[gp.X, ]
        result$VIP = VIP
    }
     
    #-- valeurs sortantes --#
    result$what = what
    result$digits = digits
    if(class(object)[1] == "pls") {
        result$method = 'pls'
	} 
	else {
        result$method = 'spls'
        result$keepX = object$keepX
        result$keepY = object$keepY
    } 
     
    class(result) = "summary"
    return(invisible(result))
}


#-------------------------- rcc -------------------------#
summary.rcc <-
function(object, 
         what = c("all", "communalities", "redundancy"), 
         cutoff = NULL, 
         digits = 4,
	 ...) 
{

    #-- initialisation des matrices --#
    ncomp = object$ncomp
    p = ncol(object$X)
    q = ncol(object$Y)
    n = nrow(object$X)
     
    bisect = object$variates$X[, 1:ncomp] + object$variates$Y[, 1:ncomp]
    cord.X = cor(object$X, bisect, use = "pairwise")
    cord.Y = cor(object$Y, bisect, use = "pairwise")
    gp.X = 1:p
    gp.Y = 1:q
     
    if (!is.null(cutoff)) {
        # choix des variables avec au moins une coordonn?e # 
        # sup?rieur au cutoff                              #
        #--------------------------------------------------#
        gp.X = vector(mode = "numeric")
        gp.Y = vector(mode = "numeric")
          
        k = 1
        for (i in 1:p) {  
            if (any(abs(cord.X[i, ]) > cutoff)) { 
                gp.X[k] = i
                k = k + 1
            }
        }
         
        k = 1
        for (i in 1:q) {  
            if (any(abs(cord.Y[i, ]) > cutoff)) { 
                gp.Y[k] = i
                k = k + 1
            }
        }
    }
     
    if (length(gp.X) == 0 || length(gp.Y) == 0) 
        stop("Cutoff value very high for the components 1:ncomp.
             No variable was selected.")
     
    result = list( )
     
    #-- communaut? --#
    if (any(what == "all") || any(what == "communalities")) {
        # communaut? inter-groupe #
        #-------------------------#
        Cm.XvsV = cor(object$X[, gp.X], object$variate$Y[, 1:ncomp], use = "pairwise")
        Cm.XvsV = t(apply(Cm.XvsV^2, 1, cumsum))
        colnames(Cm.XvsV) = paste("comp", 1:ncomp)
         
        Cm.YvsU = cor(object$Y[, gp.Y], object$variate$X[, 1:ncomp], use = "pairwise")
        Cm.YvsU = t(apply(Cm.YvsU^2, 1, cumsum))
        colnames(Cm.YvsU) = paste("comp", 1:ncomp)
         
        # communaut? intra-groupe #
        #-------------------------#
        Cm.XvsU = cor(object$X[, gp.X], object$variates$X[, 1:ncomp], use = "pairwise")
        Cm.XvsU = scale(Cm.XvsU^2, center = FALSE, scale = (object$cor[1:ncomp])^2)
        Cm.XvsU = t(apply(Cm.XvsU, 1, cumsum))
        colnames(Cm.XvsU) = paste("comp", 1:ncomp)
         
        Cm.YvsV = cor(object$Y[, gp.Y], object$variates$Y[, 1:ncomp], use = "pairwise")
        Cm.YvsV = scale(Cm.YvsV^2, center = FALSE, scale = (object$cor[1:ncomp])^2)
        Cm.YvsV = t(apply(Cm.YvsV, 1, cumsum))
        colnames(Cm.YvsV) = paste("comp", 1:ncomp)
         
        result$Cm.X = list(own = Cm.XvsU, opp = Cm.XvsV)
        result$Cm.Y = list(own = Cm.YvsV, opp = Cm.YvsU)
    }
     
    #-- redondance --#
    if (any(what == "all") || any(what == "redundancy")) {
        Rd.XvsU = cor(object$X, object$variates$X[, 1:ncomp], use = "pairwise")
        Rd.XvsU = apply(Rd.XvsU^2, 2, sum)/p
        Rd.XvsV = cor(object$X, object$variates$Y[, 1:ncomp], use = "pairwise")
        Rd.XvsV = apply(Rd.XvsV^2, 2, sum)/p
         
        Rd.YvsU = cor(object$Y, object$variates$X[, 1:ncomp], use = "pairwise")
        Rd.YvsU = apply(Rd.YvsU^2, 2, sum)/q
        Rd.YvsV = cor(object$Y, object$variates$Y[, 1:ncomp], use = "pairwise")
        Rd.YvsV = apply(Rd.YvsV^2, 2, sum)/q
         
        own = cbind(Rd.XvsU, cumsum(Rd.XvsU))
        colnames(own) = c("Proportion", "Cumulative")
        rownames(own) = paste("comp", 1:ncomp)
        opp = cbind(Rd.XvsV, cumsum(Rd.XvsV))
        colnames(opp) = c("Proportion", "Cumulative")
        rownames(opp) = paste("comp", 1:ncomp)
        result$Rd.X = list(own = own, opp = opp)
         
        own = cbind(Rd.YvsV, cumsum(Rd.YvsV))
        colnames(own) = c("Proportion", "Cumulative")
        rownames(own) = paste("comp", 1:ncomp)
        opp = cbind(Rd.YvsU, cumsum(Rd.YvsU))
        colnames(opp) = c("Proportion", "Cumulative")
        rownames(opp) = paste("comp", 1:ncomp)
        result$Rd.Y = list(own = own, opp = opp)
    }
     
    #-- corr?lations canoniques --#
    can.cor = object$cor[1:ncomp]
    names(can.cor) = paste(1:ncomp, "th", sep = "")
     
    #-- valeurs sortantes --#
    result$can.cor = can.cor
    result$ncomp = ncomp
    result$cor = can.cor
    result$cutoff = cutoff
    result$keep = list(X = colnames(object$X)[gp.X], Y = colnames(object$Y)[gp.Y])
    result$what = what
    result$digits = digits
    result$method = 'rcc'
     
    class(result) = "summary"
    return(invisible(result))
}

# from summary.prcomp, adapted for mixOmics
summary.pca <-
function (object, ...)
{
    chkDots(...)
    vars <- object$explained_variance
    importance <- rbind(`Standard deviation` = object$sdev, `Proportion of Variance` = round(vars,
    5), `Cumulative Proportion` = round(cumsum(vars), 5))
    k <- ncol(object$rotation)
    colnames(importance) <- c(colnames(object$rotation), rep("",
    length(vars) - k))
    object$importance <- importance
    class(object) <- "summary.prcomp"
    object
}

