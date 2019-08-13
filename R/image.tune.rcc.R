################################################################################
# Authors:
#   Ignacio Gonzalez,
#   Kim-Anh Le Cao,
#   Sebastien Dejean,
#
# created: 2009
# last modified: 2009
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


plot.tune.rcc <- image.tune.rcc <-

function(x, col = heat.colors, ...) 

{

    opar = par(no.readonly = TRUE)

    grid1 = x$grid1
    grid2 = x$grid2
    mat = x$mat
    nlevel = min(255, length(c(unique(mat))))
    if (nlevel / 2 - trunc(nlevel / 2) == 0 & nlevel > 6) nlevel = nlevel + 1
    col = col(nlevel)

    def.par = par(no.readonly = TRUE)
    layout(matrix(c(2, 1), ncol = 2, nrow = 1, byrow = FALSE),
    widths = c(1, 0.18))

    #-- layout 1 --#
    min.mat = min(mat)
    max.mat = max(mat)
    par(pty = "m", mai = c(1.1, 0.1, 0.95, 0.5))
    z = seq(min.mat, max.mat, length = nlevel)
    breaks = seq(min.mat, max.mat, length = nlevel + 1)
    z = matrix(z, nrow = 1)
    
    image(z, col = col, 
          zlim = c(min.mat, max.mat), oldstyle = TRUE,
          xaxt = "n", yaxt = "n", breaks = breaks)
    box()

    par(usr = c(0, 1, min.mat, max.mat))
    binwidth = (max.mat - min.mat) / nlevel
    midpoints = seq(min.mat + binwidth/2, max.mat - binwidth/2, by = binwidth)
    
    if (nlevel <= 6)
    {
        axis(4, at = midpoints, labels = round(midpoints, 3))
    }else{
        binwidth = seq(1, nlevel, by = round(nlevel / 5))
        if (binwidth[length(binwidth)] != nlevel){
            binwidth[length(binwidth) + 1] = nlevel
        }

        axis(4, at = midpoints[binwidth], labels = round(midpoints[binwidth],
        3))

    }

     

    #-- layout 2 --#
    par(mai = c(0.9, 0.85, 0.75, 0.2))
    image(grid1, grid2, mat, col = col,
    xlab = expression(lambda[1]), ylab=expression(lambda[2]),
    main = expression(CV(lambda[1], lambda[2])), axes = FALSE,
    zlim = c(min.mat, max.mat), oldstyle = TRUE)


    if (length(grid1) > 10) {
        grid1 = seq(min(grid1), max(grid1), length = 11)
    }
    if (length(grid2) > 10) {
        grid2 = seq(min(grid2), max(grid2), length = 11)
    }
    
    axis(1, at = grid1, labels = as.character(round(grid1, 4)))
    axis(2, at = grid2, labels = as.character(round(grid2, 4)))
    box()
    
    par(opar) #reset par

}

