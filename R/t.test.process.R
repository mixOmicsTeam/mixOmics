#############################################################################################################
# Author :
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
#
# created:  27-07-2017
# last modified: 27-07-2017
#
# Copyright (C) 2017
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

t.test.process = function(mat.error.rate, alpha = 0.01)
{
    # mat.error.rate has nrep rows and ncomp columns
    # we test successively whether adding a component improves the results
    
    max = ncol(mat.error.rate) #number max of components included
    pval = NULL
    opt = 1 #initialise the first optimal number of components
    for(opt in 1:max)
    {
        j=opt+1
        temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=TRUE) #t.test of "is adding X comp improves the overall results"
        if(is(temp, "try-error") || is.na(temp)) # temp can be NaN when error.keepX is constant
        {
            pval = 1
        } else {
            pval = temp
        }
        #print(opt)
        #print(j)
        #print(pval)

        while(pval> (alpha) & j<max)
        {
            j=j+1
            temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=TRUE) #t.test of "is adding X comp improves the overall results"
            if(is(temp, "try-error") || is.na(temp)) # temp can be NaN when error.keepX is constant
            {
                pval = 1
            } else {
                pval = temp
            }
            #print(opt)
            #print(j)
            #print(pval)
        }
        
        if( (pval> (alpha))) #if all pvalues were greater than alpha, then we do not increase opt and get out of the loop
        break
    }
    ncomp_opt = opt
    
    return(ncomp_opt)
}

