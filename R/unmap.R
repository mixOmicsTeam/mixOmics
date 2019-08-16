#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
#
# This function was borrowed from the mclust package and modified for mixOmics
#
# created: 2013
# last modified: 12-04-2016
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

# ---------------------------------------------------
# unmap variates.A variable for (s)plsda
# ---------------------------------------------------







#' Dummy matrix for an outcome factor
#'
#' Converts a class or group vector or factor into a matrix of indicator
#' variables.
#'
#'
#' @param classification A numeric or character vector or factor. Typically the
#' distinct entries of this vector would represent a classification of
#' observations in a data set.
#' @param groups A numeric or character vector indicating the groups from which
#' \code{classification} is drawn. If not supplied, the default is to assumed
#' to be the unique entries of classification.
#' @param noise A single numeric or character value used to indicate the value
#' of \code{groups} corresponding to noise.
#' @return An \emph{n} by \emph{K} matrix of \emph{(0,1)} indicator variables,
#' where \emph{n} is the length of samples and \emph{K} the number of classes
#' in the outcome.
#'
#' If a \code{noise} value of symbol is designated, the corresponding indicator
#' variables are relocated to the last column of the matrix.
#'
#' Note: - you can remap an unmap vector using the function \code{map} from the
#' package \pkg{mclust}. - this function should be used to unmap an outcome
#' vector as in the non-supervised methods of mixOmics. For other supervised
#' analyses such as (s)PLS-DA, (s)gccaDA this function is used internally.
#' @section References: C. Fraley and A. E. Raftery (2002). Model-based
#' clustering, discriminant analysis, and density estimation. \emph{Journal of
#' the American Statistical Association 97:611-631}.
#'
#' C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012). mclust Version
#' 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification,
#' and Density Estimation. Technical Report No. 597, Department of Statistics,
#' University of Washington.
#' @keywords cluster
#' @examples
#' \dontrun{
#' Y = unmap(nutrimouse$diet)
#' Y
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
#' # data could then used as an input in wrapper.rgcca, which is not, technically,
#' # a supervised method, see ??wrapper.rgcca
#' }
#' @export unmap
unmap = function (classification, groups = NULL, noise = NULL)
{
    n = length(classification)
    u = sort(unique(classification))
    levels =  levels(classification)### Add levels

    if (is.null(groups))
    {
        groups = u
    } else {
        if (any(match(u, groups, nomatch = 0) == 0))
        stop("groups incompatible with classification")
        miss = match(groups, u, nomatch = 0) == 0
    }

    cgroups = as.character(groups)
    if (!is.null(noise))
    {
        noiz = match(noise, groups, nomatch = 0)
        if (any(noiz == 0))
        stop("noise incompatible with classification")

        groups = c(groups[groups != noise], groups[groups == noise])
        noise = as.numeric(factor(as.character(noise), levels = unique(groups)))
    }

    groups = as.numeric(factor(cgroups, levels = unique(cgroups)))
    classification = as.numeric(factor(as.character(classification), levels = unique(cgroups)))
    k = length(groups) - length(noise)
    nam = levels(groups)

    if (!is.null(noise))
    {
        k = k + 1
        nam = nam[1:k]
        nam[k] = "noise"
    }

    z = matrix(0, n, k, dimnames = c(names(classification), nam))
    for (j in 1:k) z[classification == groups[j], j] = 1
    attr(z, "levels") = levels
    z
}

# ---------------------------------------------------
# map variable for (s)plsda
# ---------------------------------------------------







#' Classification given Probabilities
#'
#' Converts a matrix in which each row sums to \emph{1} into the nearest matrix
#' of \emph{(0,1)} indicator variables.
#'
#'
#' @param Y A matrix (for example a matrix of conditional probabilities in
#' which each row sums to 1).
#' @return A integer vector with one entry for each row of Y, in which the
#' \emph{i}-th value is the column index at which the \emph{i}-th row of
#' \code{Y} attains a maximum.
#' @section References: C. Fraley and A. E. Raftery (2002). Model-based
#' clustering, discriminant analysis, and density estimation. \emph{Journal of
#' the American Statistical Association 97:611-631}.
#'
#' C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012). mclust Version
#' 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification,
#' and Density Estimation. Technical Report No. 597, Department of Statistics,
#' University of Washington.
#' @seealso \code{\link{unmap}}
#' @examples
#'#' \dontrun{
#' Y = unmap(nutrimouse$diet)
#'
#' map(Y)
#' }
#' @export map
map = function (Y)
{
    nrowY = nrow(Y)
    cl = numeric(nrowY)
    I = 1:nrowY
    J = seq_len(Y)
    for (i in I)
    {
        cl[i] = (J[Y[i, ] == max(Y[i, ])])[1]
    }
    return(cl)
}
