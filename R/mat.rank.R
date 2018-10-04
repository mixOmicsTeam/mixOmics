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
    if (missing(tol)) 
    tol = max(dim(mat)) * max.d * .Machine$double.eps
    
    r = sum(d > tol)
     
    return(list(rank = r, tol = tol))
}

