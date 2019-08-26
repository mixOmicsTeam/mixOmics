###############################################################################
# Authors:
#  Sebastien Dejean,
#  Gonzalez, 
#  Kim-Anh Le Cao,
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
################################################################################


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








#' Matrix Rank
#' 
#' This function estimate the rank of a matrix.
#' 
#' \code{mat.rank} estimate the rank of a matrix by computing its singular
#' values \eqn{d[i]} (using \code{nipals}). The rank of the matrix can be
#' defined as the number of singular values \eqn{d[i] > 0}.
#' 
#' If \code{tol} is missing, it is given by
#' \code{tol=max(dim(mat))*max(d)*.Machine$double.eps}.
#' 
#' @param mat a numeric matrix or data frame that can contain missing values.
#' @param tol positive real, the tolerance for singular values, only those with
#' values larger than \code{tol} are considered non-zero.
#' @return The returned value is a list with components: \item{rank}{a integer
#' value, the matrix rank.} \item{tol}{the tolerance used for singular values.}
#' @author Sébastien Déjean and Ignacio González.
#' @seealso \code{\link{nipals}}
#' @keywords algebra
#' @examples
#' 
#' ## Hilbert matrix
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#' mat <- hilbert(16)
#' mat.rank(mat)
#' 
#' \dontrun{
#' ## Hilbert matrix with missing data
#' idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
#' m.na <- m <- hilbert(9)[, 1:6]
#' m.na[idx.na == 0] <- NA
#' mat.rank(m)
#' mat.rank(m.na)
#' }
#' 
#' @export mat.rank
mat.rank = function (mat, tol)
{

    if (length(dim(mat)) != 2) 
    stop("'mat' must be a numeric matrix.")
     
    mat = as.matrix(mat)
     
    if (!is.numeric(mat)) 
    stop("'mat' must be a numeric matrix.")
     
    d = nipals(mat)$eig
    max.d = d[1]
    min.d = d[length(d)]
    if (is_null(tol)) 
    tol = max(dim(mat)) * max.d * .Machine$double.eps
    
    r = sum(d > tol)
     
    return(list(rank = r, tol = tol))
}

