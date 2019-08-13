############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
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


# Copyright (C) 2009
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
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








#' Canonical Correlations Plot
#'
#' This function provides scree plot of the canonical correlations.
#'
#'
#' @param x object of class inheriting from \code{"rcc"}.
#' @param scree.type character string, (partially) matching one of
#' \code{"pointplot"} or \code{"barplot"}, determining the kind of scree plots
#' to be produced.
#' @param list() arguments to be passed to other methods. For the
#' \code{"pointplot"} type see \code{\link{points}}, for \code{"barplot"} type
#' see \code{\link{barplot}}.
#' @return none
#' @author Sébastien Déjean and Ignacio González.
#' @seealso \code{\link{points}}, \code{\link{barplot}}, \code{\link{par}}.
#' @keywords multivariate hplot
#' @examples
#'
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)
#'
#' ## 'pointplot' type scree
#' plot(nutri.res) #(default)
#'
#' \dontrun{
#' plot(nutri.res, pch = 19, cex = 1.2,
#' col = c(rep("red", 3), rep("darkblue", 18)))
#'
#' ## 'barplot' type scree
#' plot(nutri.res, scree.type = "barplot")
#'
#' plot(nutri.res, scree.type = "barplot", density = 20, col = "black")
#' }
#'
plot.rcc <-
function(x, scree.type = c("pointplot", "barplot"), ...)
{

    scree.type = match.arg(scree.type)
    if (scree.type == "pointplot") {
        plot(x$cor, xlab = "Dimension", ylim = c(0, 1),
        ylab = "Canonical correlation", ...)
    }
    else {
        barplot(x$cor, xlab = "Dimension", ylim = c(0, 1),
        ylab = "Canonical correlation", ...)
    }
}
