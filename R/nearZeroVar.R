#############################################################################################################
# Authors:
#   This function was borrowed from the package caret nzv.R with some enhancements made by
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2014
# last modified: 12-04-2016
#
# Copyright (C) 2014
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


# --------------------------------------
# nearZeroVar
# --------------------------------------







#' Identification of zero- or near-zero variance predictors
#'
#' Borrowed from the \pkg{caret} package. It is used as an internal function in
#' the PLS methods, but ca n also be used as an extermnal function, in
#' particular when the data contain a lot of zeroes values and need to be
#' prefiletered beforehand.
#'
#' This function diagnoses predictors that have one unique value (i.e. are zero
#' variance predictors) or predictors that are have both of the following
#' characteristics: they have very few unique values relative to the number of
#' samples and the ratio of the frequency of the most common value to the
#' frequency of the second most common value is large.
#'
#' For example, an example of near zero variance predictor is one that, for
#' 1000 samples, has two distinct values and 999 of them are a single value.
#'
#' To be flagged, first the frequency of the most prevalent value over the
#' second most frequent value (called the ``frequency ratio'') must be above
#' \code{freqCut}. Secondly, the ``percent of unique values,'' the number of
#' unique values divided by the total number of samples (times 100), must also
#' be below \code{uniqueCut}.
#'
#' In the above example, the frequency ratio is 999 and the unique value
#' percentage is 0.0001.
#'
#' @param x a numeric vector or matrix, or a data frame with all numeric data.
#' @param freqCut the cutoff for the ratio of the most common value to the
#' second most common value.
#' @param uniqueCut the cutoff for the percentage of distinct values out of the
#' number of total samples.
#' @return \code{nearZeroVar} returns a list that contains the following
#' components:
#'
#' \item{Position}{a vector of integers corresponding to the column positions
#' of the problematic predictors that will need to be removed.}
#' \item{Metrics}{a data frame containing the zero- or near-zero predictors
#' information with columns: \code{freqRatio}, the ratio of frequencies for the
#' most common value over the second most common value and,
#' \code{percentUnique}, the percentage of unique data points out of the total
#' number of data points.}
#' @author Max Kuhn, with speed improvements to nearZerVar by Allan Engelhardt;
#' enhancements by Florian Rohart, and speed up improvements by Benoit Gautier
#' for mixOmics
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}
#' @keywords utilities
#' @examples
#'\dontrun{
#' library(mixOmics.data)
#' nzv = nearZeroVar(diverse.16S$data.raw)
#' length(nzv$Position) # those would be removed for the default frequency cut
#'}
#' @export nearZeroVar
nearZeroVar = function (x, freqCut = 95/5, uniqueCut = 10)
{

    if (is.vector(x))
    x = matrix(x, ncol = 1)

    freqRatio = apply(x, 2, function(data)
        {
            data = na.omit(data)

            if (length(unique(data)) == length(data))
            { # No duplicate
                return(1)
            } else if (length(unique(data)) == 1) { # Same value
                return(0)
            } else {
                t = table(data)
                return(max(t, na.rm = TRUE)/max(t[-which.max(t)], na.rm = TRUE))
            }
        })

    lunique = apply(x, 2, function(data) length(unique(data[!is.na(data)])))

    percentUnique = 100 * lunique/nrow(x)
    zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))

    out = list()
    out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
    names(out$Position) = NULL
    out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
    out$Metrics = out$Metrics[out$Position, ]
    return(out)
}
