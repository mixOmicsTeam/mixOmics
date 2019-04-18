#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence in Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# This function was borrowed from the mclust package and modified for mixOmics
#
# created: 2009
# last modified: 08-07-2016
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

vip =
function(object)
{
    if (any(class(object) %in% c("mixo_plsda","mixo_splsda")))
    {
        object$Y = object$ind.mat
    } else if (any(class(object) %in% c("mixo_pls","mixo_spls"))) {
        #nothing
    } else {
        stop( " 'vip' is only implemented for the following objects: pls, plsda, spls, splsda", call.=FALSE)
    }
    #-- initialisation des matrices --#
    W = object$loadings$X
    H = object$ncomp
    q = ncol(object$Y)
    p = ncol(object$X)
    VIP = matrix(0, nrow = p, ncol = H)
    
    cor2 = cor(object$Y, object$variates$X, use = "pairwise")^2
    cor2 = as.matrix(cor2, nrow = q)
     
    VIP[, 1] = W[, 1]^2
     
    if (H > 1)
    {
        for (h in 2:H)
        {
            if (q == 1)
            {
                Rd = cor2[, 1:h] 
                VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
            } else {
                Rd = apply(cor2[, 1:h], 2, sum)
                VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
            }
        }
    }
     
    #-- valeurs sortantes --#
    VIP = sqrt(p * VIP)
    rownames(VIP) = rownames(W)
    colnames(VIP)= paste0("comp", 1:H)
     
    return(invisible(VIP))
}
