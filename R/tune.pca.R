#############################################################################################################
# Author :
#   Kim-Anh Le Cao, ARC Centre of Excellence in Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Leigh Coonan, Student, University of Quuensland, Australia
#   Fangzhou Yao, Student, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2011
# last modified: 21-04-2016
#
# Copyright (C) 2011
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







#' Tune the number of principal components in PCA
#'
#' \code{tune.pca} can be used to quickly visualise the proportion of explained
#' variance for a large number of principal components in PCA.
#'
#' The calculation is done either by a singular value decomposition of the
#' (possibly centered and scaled) data matrix, if the data is complete or by
#' using the NIPALS algorithm if there is data missing. Unlike
#' \code{\link{princomp}}, the print method for these objects prints the
#' results in a nice format and the \code{plot} method produces a bar plot of
#' the percentage of variance explaned by the principal components (PCs).
#'
#' When using NIPALS (missing values), we make the assumption that the first
#' (\code{min(ncol(X),} \code{nrow(X)}) principal components will account for
#' 100 \% of the explained variance.
#'
#' Note that \code{scale= TRUE} cannot be used if there are zero or constant
#' (for \code{center = TRUE}) variables.
#'
#' Components are omitted if their standard deviations are less than or equal
#' to \code{comp.tol} times the standard deviation of the first component. With
#' the default null setting, no components are omitted. Other settings for
#' \code{comp.tol} could be \code{comp.tol = sqrt(.Machine$double.eps)}, which
#' would omit essentially constant components, or \code{comp.tol = 0}.
#'
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#'
#' @param X a numeric matrix (or data frame) which provides the data for the
#' principal components analysis. It can contain missing values.
#' @param ncomp integer, the number of components to initially analyse in
#' \code{tune.pca} to choose a final \code{ncomp} for \code{pca}. If
#' \code{NULL}, function sets \code{ncomp = min(nrow(X), ncol(X))}
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
#' @param max.iter integer, the maximum number of iterations for the NIPALS
#' algorithm.
#' @param tol a positive real, the tolerance used for the NIPALS algorithm.
#' @param logratio one of ('none','CLR','ILR'). Default to 'none'
#' @param V Matrix used in the logratio transformation id provided.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements).
#' @return \code{tune.pca} returns a list with class \code{"tune.pca"}
#' containing the following components: \item{sdev}{the square root of the
#' eigenvalues of the covariance/correlation matrix, though the calculation is
#' actually done with the singular values of the data matrix).}
#' \item{explained_variance}{the proportion of explained variance accounted for
#' by each principal component is calculated using the eigenvalues}
#' \item{cum.var}{the cumulative proportion of explained variance accounted for
#' by the sequential accumulation of principal components is calculated using
#' the sum of the proportion of explained variance}
#' @author Ignacio González and Leigh Coonan
#' @seealso \code{\link{nipals}}, \code{\link{biplot}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}} and http://www.mixOmics.org
#' for more details.
#' @keywords algebra
#' @examples
#' \dontrun{
#' library(mixOmics.data)
#' tune <- tune.pca(liver.toxicity$gene, center = TRUE, scale = TRUE)
#' tune
#'}
#' @export tune.pca
tune.pca =
function(X,
ncomp = NULL,
center = TRUE, 	# sets the mean of the data to zero, ensures that the first PC describes the direction of the maximum variance
scale = FALSE, 	# variance is unit accross different units
max.iter = 500,
tol = 1e-09,
logratio = 'none',# one of ('none','CLR','ILR')
V = NULL,
multilevel = NULL)
{


    result = pca(X = X, ncomp = ncomp,
    center = center, scale = scale,
    max.iter = max.iter, tol = tol,
    logratio = logratio, V = V,
    multilevel = multilevel)

    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE

    #  list eigenvalues, prop. of explained varience and cumulative proportion of explained variance
    prop.var = result$explained_variance
    cum.var = result$cum.var

    ind.show = min(10, ncomp)

    print(result)

    # Plot the principal components and explained variance
    # note: if NA values, we have an estimation of the variance using NIPALS
    if(!na.X)
    {
        ylab = "Proportion of Explained Variance"
    } else{
        ylab = "Estimated Proportion of Explained Variance"
    }
    barplot(prop.var[1:result$ncomp], names.arg = 1:result$ncomp, xlab = "Principal Components",
    ylab = ylab)

    result$call = match.call()

    class(result) = "tune.pca"
    return(invisible(result))
}
