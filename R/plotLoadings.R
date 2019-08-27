#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 15-04-2016
# last modified: 24-05-2016
#
# Copyright (C) 2016
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








#' Plot of Loading vectors
#'
#' This function provides a horizontal bar plot to visualise loading vectors.
#' For discriminant analysis, it provides visualisation of highest or lowest
#' mean/median value of the variables with color code corresponding to the
#' outcome of interest.
#'
#'
#' The contribution of each variable for each component (depending on the
#' object) is represented in a barplot where each bar length corresponds to the
#' loading weight (importance) of the feature. The loading weight can be
#' positive or negative.
#'
#' For discriminant analysis, the color corresponds to the group in which the
#' feature is most 'abundant'.  Note that this type of graphical output is
#' particularly insightful for count microbial data - in that latter case using
#' the \code{method = 'median'} is advised. Note also that if the parameter
#' \code{contrib} is not provided, plots are white.
#'
#' For MINT analysis, \code{study="global"} plots the global loadings while
#' partial loadings are plotted when \code{study} is a level of
#' \code{object$study}. Since variable selection in MINT is performed at the
#' global level, only the selected variables are plotted for the partial
#' loadings even if the partial loadings are not sparse. See references.
#' Importantly for multi plots, the legend accounts for one subplot in the
#' layout design.
#'
#' @aliases plotLoadings plotLoadings.mixo_pls plotLoadings.mixo_spls
#' plotLoadings.rcc plotLoadings.pca plotLoadings.sgcca plotLoadings.rgcca
#' plotLoadings.block.pls plotLoadings.block.spls plotLoadings.mixo_plsda
#' plotLoadings.mixo_splsda plotLoadings.sgccda plotLoadings.block.plsda
#' plotLoadings.block.splsda plotLoadings.mint.pls plotLoadings.mint.spls
#' plotLoadings.mint.plsda plotLoadings.mint.splsda
#' @param object object
#' @param contrib a character set to 'max' or 'min' indicating if the color of
#' the bar should correspond to the group with the maximal or minimal
#' expression levels / abundance.
#' @param method a character set to 'mean' or 'median' indicating the criterion
#' to assess the contribution. We recommend using median in the case of count
#' or skewed data.
#' @param study Indicates which study are to be plotted.  A character vector
#' containing some levels of \code{object$study}, "all.partial" to plot all
#' studies or "global" is expected.
#' @param block A single value indicating which block to consider in a
#' \code{sgccda} object.
#' @param comp integer value indicating the component of interest from the
#' object.
#' @param col color used in the barplot, only for object from non Discriminant
#' analysis
#' @param plot Boolean indicating of the plot should be output. If set to FALSE
#' the user can extract the contribution matrix, see example. Default value is
#' TRUE.
#' @param show.ties Boolean. If TRUE then tie groups appear in the color set by
#' \code{col.ties}, which will appear in the legend. Ties can happen when
#' dealing with count data type. By default set to TRUE.
#' @param col.ties Color corresponding to ties, only used if
#' \code{show.ties=TRUE} and ties are present.
#' @param ndisplay integer indicating how many of the most important variables
#' are to be plotted (ranked by decreasing weights in each PLS-component).
#' Useful to lighten a graph.
#' @param size.name A numerical value giving the amount by which plotting the
#' variable name text should be magnified or reduced relative to the default.
#' @param size.legend A numerical value giving the amount by which plotting the
#' legend text should be magnified or reduced relative to the default.
#' @param name.var A character vector indicating the names of the variables.
#' The names of the vector should match the names of the input data, see
#' example.
#' @param name.var.complete Boolean. If \code{name.var} is supplied with some
#' empty names, \code{name.var.complete} allows you to use the initial variable
#' names to complete the graph (from colnames(X)). Defaut to FALSE.
#' @param title A set of characters to indicate the title of the plot. Default
#' value is NULL.
#' @param subtitle subtitle for each plot, only used when several \code{block}
#' or \code{study} are plotted.
#' @param size.title size of the title
#' @param size.subtitle size of the subtitle
#' @param legend Boolean indicating if the legend indicating the group outcomes
#' should be added to the plot. Default value is TRUE.
#' @param legend.color A color vector of length the number of group outcomes.
#' See examples.
#' @param legend.title A set of characters to indicate the title of the legend.
#' Default value is NULL.
#' @param layout Vector of two values (rows,cols) that indicates the layout of
#' the plot. If \code{layout} is provided, the remaining empty subplots are
#' still active
#' @param border Argument from \code{\link{barplot}}: indicates whether to draw
#' a border on the barplot.
#' @param xlim Argument from \code{\link{barplot}}: limit of the x-axis. When
#' plotting several \code{block}, a matrix is expected where each row is the
#' \code{xlim} used for each of the blocks.
#' @param \dots not used.
#' @return none
#' @author Florian Rohart, Kim-Anh Lê Cao, Benoit Gautier
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{mint.pls}}, \code{\link{mint.spls}},
#' \code{\link{mint.plsda}}, \code{\link{mint.splsda}},
#' \code{\link{block.pls}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}}, \code{\link{block.splsda}},
#' \code{\link{mint.block.pls}}, \code{\link{mint.block.spls}},
#' \code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}}
#' @references Rohart F. et al (2016, submitted). MINT: A multivariate
#' integrative approach to identify a reproducible biomarker signature across
#' multiple experiments and platforms.
#'
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group
#' PLS Regression: Application to Epidemiology. In New Perspectives in Partial
#' Least Squares and Related Methods, pages 243-255. Springer.
#'
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO - multi omics integration for biomarker discovery.
#'
#' Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
#' Sparse canonical methods for biological data integration: application to a
#' cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#'
#' Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris:
#' Editions Technic.
#'
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' @keywords multivariate
#' @examples
#'
#' \dontrun{
#' ## object of class 'spls'
#' # --------------------------
#' library(mixOmicsData)
#' X = liver.toxicity$gene
#' Y = liver.toxicity$clinic
#'
#' toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50),
#' keepY = c(10, 10))
#'
#' plotLoadings(toxicity.spls)
#'
#' # with xlim
#' xlim = matrix(c(-0.1,0.3, -0.4,0.6), nrow = 2, byrow = TRUE)
#' plotLoadings(toxicity.spls, xlim = xlim)
#'
#'
#' ## object of class 'splsda'
#' # --------------------------
#' X = as.matrix(liver.toxicity$gene)
#' Y = as.factor(liver.toxicity$treatment[, 4])
#'
#' splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))
#'
#' # contribution on comp 1, based on the median.
#' # Colors indicate the group in which the median expression is maximal
#' plotLoadings(splsda.liver, comp = 1, method = 'median')
#' plotLoadings(splsda.liver, comp = 1, method = 'median', contrib = "max")
#'
#' # contribution on comp 2, based on median.
#' #Colors indicate the group in which the median expression is maximal
#' plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = "max")
#'
#' # contribution on comp 2, based on median.
#' # Colors indicate the group in which the median expression is minimal
#' plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = 'min')
#'
#' # changing the name to gene names
#' # if the user input a name.var but names(name.var) is NULL,
#' # then a warning will be output and assign names of name.var to colnames(X)
#' # this is to make sure we can match the name of the selected variables to the contribution plot.
#' name.var = liver.toxicity$gene.ID[, 'geneBank']
#' length(name.var)
#' plotLoadings(splsda.liver, comp = 2, method = 'median', name.var = name.var,
#' title = "Liver data", contrib = "max")
#'
#' # if names are provided: ok, even when NAs
#' name.var = liver.toxicity$gene.ID[, 'geneBank']
#' names(name.var) = rownames(liver.toxicity$gene.ID)
#' plotLoadings(splsda.liver, comp = 2, method = 'median',
#' name.var = name.var, size.name = 0.5, contrib = "max")
#'
#' #missing names of some genes? complete with the original names
#' plotLoadings(splsda.liver, comp = 2, method = 'median',
#' name.var = name.var, size.name = 0.5,complete.name.var=TRUE, contrib = "max")
#'
#' # look at the contribution (median) for each variable
#' plot.contrib = plotLoadings(splsda.liver, comp = 2, method = 'median', plot = FALSE,
#' contrib = "max")
#' head(plot.contrib$contrib)
#' # change the title of the legend and title name
#' plotLoadings(splsda.liver, comp = 2, method = 'median', legend.title = 'Time',
#' title = 'Contribution plot', contrib = "max")
#'
#' # no legend
#' plotLoadings(splsda.liver, comp = 2, method = 'median', legend = FALSE, contrib = "max")
#'
#' # change the color of the legend
#' plotLoadings(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4), contrib = "max")
#'
#'
#'
#' # object 'splsda multilevel'
#' # -----------------
#'
#' X = vac18$genes
#' Y = vac18$stimulation
#' # sample indicates the repeated measurements
#' sample = vac18$sample
#' stimul = vac18$stimulation
#'
#' # multilevel sPLS-DA model
#' res.1level = splsda(X, Y = stimul, ncomp = 3, multilevel = sample,
#' keepX = c(30, 137, 123))
#'
#'
#' name.var = vac18$tab.prob.gene[, 'Gene']
#' names(name.var) = colnames(X)
#'
#' plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
#' name.var = name.var, size.name = 0.2, contrib = "max")
#'
#' # too many transcripts? only output the top ones
#' plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
#' name.var = name.var, size.name = 0.5, ndisplay = 60, contrib = "max")
#'
#'
#'
#' # object 'plsda'
#' # ----------------
#'
#' # breast tumors
#' # ---
#' X = breast.tumors$gene.exp
#' Y = breast.tumors$sample$treatment
#'
#' plsda.breast = plsda(X, Y, ncomp = 2)
#'
#' name.var = as.character(breast.tumors$genes$name)
#' names(name.var) = colnames(X)
#'
#' # with gene IDs, showing the top 60
#' plotLoadings(plsda.breast, contrib = 'max', comp = 1, method = 'median',
#' ndisplay = 60,
#' name.var = name.var,
#' size.name = 0.6,
#' legend.color = color.mixo(1:2))
#'
#'
#' # liver toxicity
#' # ---
#'
#' X = liver.toxicity$gene
#' Y = liver.toxicity$treatment[, 4]
#'
#' plsda.liver = plsda(X, Y, ncomp = 2)
#' plotIndiv(plsda.liver, ind.names = Y, ellipse = TRUE)
#'
#'
#' name.var = liver.toxicity$gene.ID[, 'geneBank']
#' names(name.var) = rownames(liver.toxicity$gene.ID)
#'
#' plotLoadings(plsda.liver, contrib = 'max', comp = 1, method = 'median', ndisplay = 100,
#' name.var = name.var, size.name = 0.4,
#' legend.color = color.mixo(1:4))
#'
#'
#' # object 'sgccda'
#' # ----------------
#'
#' Y = nutrimouse$diet
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
#' design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#'
#' nutrimouse.sgccda = wrapper.sgccda(X = data,
#' Y = Y,
#' design = design,
#' keepX = list(gene = c(10,10), lipid = c(15,15)),
#' ncomp = 2,
#' scheme = "centroid")
#'
#' plotLoadings(nutrimouse.sgccda,block=2)
#' plotLoadings(nutrimouse.sgccda,block="gene")
#'
#'
#'
#' # object 'mint.splsda'
#' # ----------------
#' data = stemcells$gene
#' type.id = stemcells$celltype
#' exp = stemcells$study
#'
#' res = mint.splsda(X = data, Y = type.id, ncomp = 3, keepX = c(10,5,15), study = exp)
#'
#' plotLoadings(res)
#' plotLoadings(res, contrib = "max")
#' plotLoadings(res, contrib = "min", study = 1:4,comp=2)
#'
#' # combining different plots by setting a layout of 2 rows and 4columns.
#' # Note that the legend accounts for a subplot so 4columns instead of 2.
#' plotLoadings(res,contrib="min",study=c(1,2,3),comp=2, layout = c(2,4))
#' plotLoadings(res,contrib="min",study="global",comp=2)
#'
#' }
#'
#' @export plotLoadings
plotLoadings  =
function(object, ...) UseMethod("plotLoadings")


# --------------------------------------------------------------------------------------
# Internal helpers functions to run "plotLoadings" functions
# --------------------------------------------------------------------------------------


check.input.plotLoadings = function(object, block, study, subtitle, size.name, size.legend, title, col, contrib, name.var, xlim)
{

    if (is.null(object$loadings))
    stop("'plotLoadings' should be used on object for which object$loadings is present.")

    # block
    # --
    if (isNULL(block))
    {
        if (!is(object, "DA"))
        {
            block = object$names$blocks
        } else  if (is(object, c("mixo_plsda", "mixo_splsda"))) {
            block = "X"
        } else {
            if (!is.null(object$indY))
            {
                block = object$names$blocks[-object$indY]
            } else {
                block = object$names$blocks
            }
        }
    }

    if (is(object, c("mixo_plsda", "mixo_splsda")) & (!all(block %in% c(1,"X")) | length(block) > 1 ))
    stop("'block' can only be 'X' or '1' for plsda and splsda object")

    if (is(object, c("mixo_plsda", "mixo_splsda","pca")))
    {
        object$indY = 2
    } else if (is(object, c("mixo_pls", "mixo_spls"))) {
        object$indY = 3 # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2
    }

    if(!is(object, "DA"))
    object$indY = length(object$names$blocks)+1  # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2

    if(is.numeric(block))
    {
        if(any(block>length(object$names$blocks[-object$indY])))
        stop("'block' needs to be lower than the number of blocks in the fitted model, which is ",length(object$names$blocks)-1)

    }else if(is.character(block) & any(is.na(match(block,object$names$blocks[-object$indY])))) {
        stop("Incorrect value for 'block', 'block' should be among the blocks used in your object: ", paste(object$names$blocks[-object$indY],collapse=", "), call. = FALSE)
    }


    if (!isNULL(subtitle))
    {
        if (length(subtitle)!=length(block))
        stop("'subtitle' indicates the subtitle of the plot for each block and it needs to be the same length as 'block'.")
    }

    if (!isNULL(study))
    {
    #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
        if (any(!study%in%c(levels(object$study), "global")))
        stop("'study' must be one of 'object$study' or 'all'.")

        if (length(study)!=length(unique(study)))
        stop("Duplicate in 'study' not allowed")
    }

    # cex
    # --
    if (size.name <= 0)
    size.name = 0.7

    if (!isNULL(size.legend))
    {
        if(size.legend <= 0)
        size.legend = 0.8
    } else {
        size.legend = NULL
    }

    # contrib
    # --
    if (!isNULL(contrib))
    {
        if(length(contrib) > 1 | !all(contrib %in% c("min", "max")))
        stop("'contrib' must be either 'min' or 'max'")

    }

    # xlim
    #---
    if (!isNULL(xlim))
    {
        # check xlim, has to be a matrix with number of rows=number of blocks, or a vector of two values
        if(length(block) == 1 & !is.null(xlim))
        {
            if(length(xlim) !=2)
            stop("'xlim' must be a vector of length 2")

            xlim = matrix(xlim, nrow = 1)
        }

        if(length(block)>1 & !is.null(xlim))
        {
            if(is.matrix(xlim) && ( !nrow(xlim) %in%c(1, length(block))  | ncol(xlim) != 2 ))
            stop("'xlim' must be a matrix with ",length(block)," rows (length(block)) and 2 columns")

            if(is.vector(xlim))
            {
                if(length(xlim) !=2)
                stop("'xlim' must be a matrix with ",length(block)," rows (length(block)) and 2 columns")

                xlim = matrix(xlim, nrow = 1)
            }

            if(nrow(xlim) != length(block)) # we complete xlim to have one xlim per block
            xlim = matrix(rep(xlim, length(block)), nrow = length(block), byrow=TRUE)
        }

    } else {
        xlim = NULL
    }

    #names.var
    #-----
    if(!is.null(name.var))
    {
        if (length(block) >1 && length(block) != length(name.var))
        stop("'names' has to be a list of length the number of block to plot: ", length(block))

        if (length(block) > 1)
        {
            for (block_i in block)
            {
                if(length(name.var[[block_i]])!= nrow(object$loadings[[block_i]]))
                stop("For block '", block_i,"', 'name.var' should be a vector of length ", nrow(object$loadings[[block_i]]))
            }
        } else {
            if(length(name.var)!= nrow(object$loadings[[block]]))
            stop("For block '", block,"', 'name.var' should be a vector of length ", nrow(object$loadings[[block]]))

        }
    }
    #title
    #-----
    if (!is.null(title) & !is.character(title))
    warning('title needs to be of type character')

    #col
    #-----
    if (!is.null(col) & (length(col) !=  1))
    {
        warning('col must be the of length 1, by default set to default colors')
        col = color.mixo(1)  # by default set to the colors in color.mixo (10 colors)
    }
    if (is.null(col))
    col = color.mixo(1) # by default set to the colors in color.mixo (10 colors)


    return(list(col = col, size.name = size.name, size.legend = size.legend, block = block, xlim = xlim))
}

layout.plotLoadings = function(layout, plot, legend, block)
{
    # layout
    # --
    if(plot == TRUE)
    {
        opar = par(no.readonly = TRUE)
        reset.mfrow = FALSE # if set to TRUE, the algorithm ends up with  par(mfrow=reset.mfrow)
        nResp = length(block) + length(block) * legend  #number of blocks *2 if legend is plotted

        if (is.null(layout))
        {
            # check if there are enough plots in mfrow
            omfrow = par("mfrow")
            available.plots = prod(omfrow)
            if (available.plots<nResp) # if not enough plots available, we create our new plot
            {

                if (legend)
                {
                    nRows = min(c(2, ceiling(nResp / 4)))
                    nCols = min(c(4, nResp))
                    layout(matrix(1 : (nCols * nRows), nRows, nCols, byrow=TRUE),rep(c(0.7,0.7 -0.4*legend),nCols/(1+legend)))
                } else {
                    nRows = min(c(3, ceiling(nResp/3)))
                    nCols = min(c(3, ceiling(nResp / nRows)))

                    layout(matrix(1 : (nCols * nRows), nRows, nCols, byrow=TRUE))

                }
                if (nRows * nCols < nResp)
                devAskNewPage(TRUE)

                reset.mfrow=TRUE # we changed mfrow to suits our needs, so we reset it at the end
            }
        } else {
            if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
            stop("'layout' must be a numeric vector of length 2.")

            nRows = layout[1]
            nCols = layout[2]
            par(mfrow = layout)

            if (nRows * nCols < nResp)
            devAskNewPage(TRUE)
        }

    } else {
        reset.mfrow = FALSE
        opar = NULL
    }

    return(list(reset.mfrow = reset.mfrow, opar = opar))

}


get.loadings.ndisplay = function(object, comp, block, name.var, name.var.complete, ndisplay)
{
    ##selectvar
    selected.var = selectVar(object, comp = comp, block = block) # gives name and values of the blocks in 'block'
    name.selected.var = selected.var[[1]]$name
    value.selected.var = selected.var[[1]]$value

    # ndisplay
    # ------
    # if null set by default to all variables from selectVar
    if (is.null(ndisplay))
    {
        ndisplay.temp = length(name.selected.var)
    } else if (ndisplay > length(name.selected.var)) {
        message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var), " for block ", block)
        ndisplay.temp = length(name.selected.var)
    } else {
        ndisplay.temp = ndisplay
    }

    name.selected.var = name.selected.var[1:ndisplay.temp]
    value.selected.var = value.selected.var[1:ndisplay.temp,]

    #comp
    # ----
    if (is(object, c("mixo_pls","mixo_spls", "rcc")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
    {
        ncomp = object$ncomp
        object$X = list(X = object$X, Y = object$Y) # so that the data is in object$X, either it's a pls or block approach
    } else {
        ncomp = object$ncomp[block]
    }

    if (any(max(comp) > ncomp))
    stop(paste("Argument 'comp' should be less or equal to ", ncomp))

    names.block = as.character(names(selected.var)[1]) #it should be one block and ncomp, so we take the first one

    X = object$X[names.block][[1]]

    #name.var
    ind.match = match(name.selected.var, colnames(X)) # look at the position of the selected variables in the original data X
    if(!is.null(name.var))
    {
        if(length(name.var)!= ncol(X))
        stop("For block '", names.block,"', 'name.var' should be a vector of length ", ncol(X))

        colnames.X = as.character(name.var[ind.match]) # get the
    }else{
        colnames.X = as.character(colnames(X))[ind.match]
    }
    X = X[, name.selected.var, drop = FALSE] #reduce the problem to ndisplay

    #completing colnames.X by the original names of the variables when missing
    if (name.var.complete == TRUE)
    {
        ind = which(colnames.X == "")
        if (length(ind) > 0)
        colnames.X[ind] = colnames(X)[ind]
    }


    return(list(X = X, names.block = names.block, colnames.X = colnames.X, name.selected.var = name.selected.var, value.selected.var = value.selected.var))
}



get.contrib.df = function(Y, X, method, contrib, value.selected.var, colnames.X, name.selected.var, legend.color, col.ties)
{
    # Start: Initialisation
    which.comp = method.group = list()
    which.contrib = data.frame(matrix(FALSE, ncol = nlevels(Y) + 2, nrow = length(colnames.X),
    dimnames = list(name.selected.var, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
    # End: Initialisation

    # calculate the max.method per group for each variable, and identifies which group has the max max.method
    for(k in 1:ncol(X))
    {
        method.group[[k]] = tapply(X[, k], Y, method, na.rm=TRUE) #method is either mean or median
        # determine which group has the highest mean/median
        which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == get(contrib)((method.group[[k]])) # contrib is either min or max
    }

    # we also add an output column indicating the group that is max
    # if ties, we set the color to white
    which.contrib$color = apply(which.contrib, 1, function(x)
    {
        if (length(which(x)) > 1)
        {
            return(col.ties)
        } else { # otherwise we use legend color provided
            return(legend.color[1 : nlevels(Y)][which(x)])
        }
    })

    which.contrib$GroupContrib = apply(which.contrib[, 1:(nlevels(Y))], 1, function(x)
    {
        if (length(which(x)) > 1)
        {
            return("tie")
        } else {
            return(levels(Y)[which(x)])
        }
    })

    method.group = do.call(rbind, method.group)
    df = data.frame(method.group, which.contrib, importance = value.selected.var)
    return(df)
}
