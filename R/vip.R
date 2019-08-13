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







#' Variable Importance in the Projection (VIP)
#'
#' The function \code{vip} computes the influence on the \eqn{Y}-responses of
#' every predictor \eqn{X} in the model.
#'
#' Variable importance in projection (VIP) coefficients reflect the relative
#' importance of each \eqn{X} variable for each \eqn{X} variate in the
#' prediction model. VIP coefficients thus represent the importance of each
#' \eqn{X} variable in fitting both the \eqn{X}- and \eqn{Y}-variates, since
#' the \eqn{Y}-variates are predicted from the \eqn{X}-variates.
#'
#' VIP allows to classify the \eqn{X}-variables according to their explanatory
#' power of \eqn{Y}. Predictors with large VIP, larger than 1, are the most
#' relevant for explaining \eqn{Y}.
#'
#' @param object object of class inheriting from \code{"pls"}, \code{"plsda"},
#' \code{"spls"} or \code{"splsda"}.
#' @return \code{vip} produces a matrix of VIP coefficients for each \eqn{X}
#' variable (rows) on each variate component (columns).
#' @author Sébastien Déjean and Ignacio Gonz\`alez.
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{summary}}.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' @keywords regression multivariate
#' @examples
#'
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y)
#'
#' linn.vip <- vip(linn.pls)
#'
#' barplot(linn.vip,
#' beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan"),
#' ylim = c(0, 1.7), legend = rownames(linn.vip),
#' main = "Variable Importance in the Projection", font.main = 4)
#'
#' @export vip
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
