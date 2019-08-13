#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
#
# created: 2015
# last modified: 17-03-2016
#
# Copyright (C) 2015
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


# ========================================================================================================
# selectVar: output the variables that were selected on 'comp', for 'block'
# ========================================================================================================

# object: a pls, spls, block, mint, rcc, sgcca  object
# comp: to display the variables selected on dimension 'comp'
# block: display the selected variables on the data 'block', from object$names$blocks







#' Output of selected variables
#'
#' This function outputs the selected variables on each component for the
#' sparse versions of the approaches (was also generalised to the non sparse
#' versions for our internal functions).
#'
#' \code{selectVar} provides the variables selected on a given component. \
#' \describe{ \item{list("name")}{outputs the name of the selected variables
#' (provided that the input data have colnames) ranked in decreasing order of
#' importance.} \item{list("value")}{outputs the loading value for each
#' selected variable, the loadings are ranked according to their absolute
#' value.} } These functions are only implemented for the sparse versions.
#'
#' @aliases selectVar selectVar.mixo_pls selectVar.mixo_spls selectVar.pca
#' selectVar.sgcca selectVar.rgcca select.var
#' @param object object of class inherited from \code{"pls"}, \code{"spls"},
#' \code{"plsda"},\code{"splsda"}, \code{"pca"}, \code{"spca"}, \code{"sipca"}.
#' @param comp integer value indicating the component of interest.
#' @param block for an object of class \code{"sgcca"}, the block data sets can
#' be specified as an input vector, for example \code{c(1,2)} for the first two
#' blocks. Default to NULL (all block data sets)
#' @param list() other arguments.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart.
#' @examples
#'
#'
#' X = liver.toxicity$gene
#' Y = liver.toxicity$clinic
#'
#' # example with sPCA
#' # ------------------
#' liver.spca <- spca(X, ncomp = 1, keepX = 10)
#' selectVar(liver.spca, comp = 1)$name
#' selectVar(liver.spca, comp = 1)$value
#'
#' \dontrun{
#' #example with sIPCA
#' # -----------------
#'
#' liver.sipca <- sipca(X, ncomp = 3, keepX = rep(10, 3))
#' selectVar(liver.sipca, comp = 1)
#'
#' # example with sPLS
#' # -----------------
#'
#' liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
#' selectVar(liver.spls, comp = 2)
#'
#' # example with sPLS-DA
#' X = srbct$gene
#' Y = srbct$class
#'
#' srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
#' select = selectVar(srbct.splsda, comp = 2)
#' select
#' # this is a very specific case where a data set has no rownames.
#' srbct$gene.name[substr(select$select, 2,5),]
#'
#'
#' # example with sGCCA
#' # -----------------
#' # ! need to unmap the Y factor
#' Y = unmap(nutrimouse$diet)
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid,Y)
#' # in this design, gene expression and lipids are connected to the diet factor
#' # and gene expression and lipids are also connected
#' design = matrix(c(0,1,1,
#' 1,0,1,
#' 1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#' #note: the penalty parameters need to be tuned
#' wrap.result.sgcca = wrapper.sgcca(X = data, design = design, penalty = c(.3,.3, 1),
#' ncomp = 2,
#' scheme = "horst")
#'
#' #variables selected and loadings values on component 1 for the two blocs
#' selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))
#'
#' #variables selected on component 1 for each block
#' selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$'gene'$name
#' selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$'lipid'$name
#'
#' #variables selected on component 2 for each block
#' selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$'gene'$name
#' selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$'lipid'$name
#'
#' # loading value of the variables selected on the first block
#' selectVar(wrap.result.sgcca, comp = 1, block = 1)$'gene'$value
#' }
#'
#' @export selectVar
selectVar <-
function(...) UseMethod("selectVar")


get.name.and.value=function(x,comp)
{
    if(length(x[,comp,drop=FALSE]) > 1)
    {
        name.var = names(sort(abs(x[,comp]), decreasing = TRUE)[1:sum(x[,comp]!=0)])
    } else {
        name.var = rownames(x) # when only one number, sort loses the name of the variable
    }
    value.var=x[name.var,comp]
    return(list(name = name.var, value = data.frame(value.var)))
}


# ------------------ for all object  --------------------
selectVar.mixo_spls  = selectVar.mixo_pls  =
selectVar.sgcca = selectVar.rgcca =
selectVar.pca =
function(object, comp =1, block=NULL, ...)
{

    # check arguments
    # -----------------
    if (length(comp) > 1)
        stop("Expecting one single value for 'comp'")

    if (is.null(block))
    {
        if (any(comp > object$ncomp))
        stop("'comp' is greater than the number of components in the fitted model")
        null.block=TRUE
        block=1:length(object$loadings)

    }else{
        if (any(class(object)%in%c("pca")))
        object$names$blocks="X"

        if (is.numeric(block))
        {
            if (any(block>length(object$names$blocks)))
            stop("'block' needs to be lower than the number of blocks in the fitted model, which is length(object$names$blocks)")

        }else if (is.character(block) & sum(!is.na(match(block,object$names$blocks)))==0) {
            stop("No entry of 'block'  match object$names$blocks")

        }else if (is.character(block) & sum(is.na(match(block,object$names$blocks)))>0) {
            warning("At least one entry of 'block' does not match object$names$blocks")
        }

        if (length(object$ncomp)>1)
        {
            if (any(comp > object$ncomp[block]))
            stop("'comp' is greater than the number of components in the fitted model for the block you specified. See object$ncomp")

        }else{
            if (any(comp > object$ncomp))
            stop("'comp' is greater than the number of components in the fitted model")
        }

        null.block=FALSE
    }

    # main function: get the names and values of the non zero loadings
    # -----------------
    out = lapply(object$loadings[block],get.name.and.value,comp=comp)


    # outputs
    # ----------
    #if all blocks are considered by default (null.block=TRUE) and it's a DA analysis, then we don't show Y
    if (null.block)
    {
        if (any(class(object)%in%c("block.plsda","block.splsda")))# the position of Y is in indY
        {
            out=out[-object$indY] #remove Y
        }else if (any(class(object)%in%c("mint.plsda","mint.splsda","mixo_plsda","mixo_splsda"))) {
            # Y is always in second position
            out=out[[1]]
        }else if (any(class(object)%in%c("pca"))) { #keep the result as a list
            out=out[[1]]
        }

    } else {
        if (length(grep("pca",class(object)))>0)
        out=out
    }

    #we add comp as an output
    out$comp=comp

    return(out)
}

