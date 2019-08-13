################################################################################
# Authors:
#   Amrit Singh,
#   Michael Vacher,
#   Florian Rohart,
#   Kim-Anh Le Cao,
#
# created: 2015
# last modified: 24-08-2016
#
# Copyright (C) 2015
#
# This program is free software  you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation  either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY  without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program  if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
################################################################################








#' circosPlot for DIABLO
#'
#' Displays variable correlation among different blocks
#'
#' \code{circosPlot} function depicts correlations of variables selected with
#' \code{block.splsda} among different blocks, using a generalisation of the
#' method presented in González et al 2012. If \code{ncomp} is specified, then
#' only the variables selected on that component are displayed.
#'
#' @param object An object of class inheriting from \code{"block.splsda"}.
#' @param comp Numeric vector indicating which component to plot. Default to
#' all
#' @param cutoff Only shows links with a correlation higher than \code{cutoff}
#' @param color.Y a character vector of colors to be used for the levels of the
#' outcome
#' @param color.blocks a character vector of colors to be used for the blocks
#' @param color.cor a character vector of two colors. First one is for the
#' negative correlation, second one is for the positive correlation
#' @param var.names Optional parameter. A list of length the number of blocks
#' in \code{object$X}, containing the names of the variables of each block. If
#' \code{NULL}, the colnames of the data matrix are used.
#' @param showIntraLinks if TRUE, shows the correlation higher than the
#' threshold inside each block.
#' @param line if TRUE, shows the overall expression of the selected variables.
#' see examples.
#' @param size.legend size of the legend
#' @param ncol.legend number of columns for the legend
#' @param size.variables size of the variable labels
#' @param size.labels size of the block labels
#' @param legend boolean. Whether the legend should be added. Default is TRUE.
#' @return If saved in an object, the circos plot will output the similarity
#' matrix and the names of the variables displayed on the plot (see
#' \code{attributes(object)}).
#' @author Michael Vacher, Amrit Singh, Florian Rohart, Kim-Anh Lê Cao
#' @seealso \code{\link{block.splsda}}, references and
#' http://www.mixOmics.org/mixDIABLO for more details.
#' @references Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt
#' S. and Lê Cao K.A. (2016). DIABLO: multi omics integration for biomarker
#' discovery. BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611}
#'
#' mixOmics article:
#'
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#'
#' González I., Lê Cao K.A., Davis M.J., Déjean S. (2012). Visualising
#' associations between paired 'omics' data sets. \emph{BioData Mining};
#' \bold{5}(1).
#' @keywords regression multivariate
#' @examples
#'
#' Y = nutrimouse$diet
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
#' design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#'
#'
#' nutrimouse.sgccda <- wrapper.sgccda(X=data,
#' Y = Y,
#' design = design,
#' keepX = list(gene=c(10,10), lipid=c(15,15)),
#' ncomp = 2,
#' scheme = "horst")
#'
#' circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1)
#'
#' \dontrun{
#' circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1,
#' color.Y = 1:5, color.blocks = c("green","brown"), color.cor = c("magenta", "purple"))
#'
#' par(mfrow=c(2,2))
#' circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2,
#' size.legend = 1.1)
#' circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2,
#' size.legend = 1.1, showIntraLinks = TRUE)
#' circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 1,
#' size.legend = 1.1, showIntraLinks = TRUE)
#' circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2,
#' size.legend = 1.1, showIntraLinks = TRUE, line = FALSE, size.variables = 0.5)
#' }
#'
#' @importFrom reshape2 dcast
#' @importFrom grDevices col2rgb
#' @export circosPlot

circosPlot = function(object,
comp = 1 : min(object$ncomp),
cutoff,
color.Y,
color.blocks,
color.cor,
var.names = NULL,
showIntraLinks = FALSE,
line=TRUE,
size.legend=0.8,
ncol.legend=1,
size.variables = 0.25,
size.labels=1,
legend = TRUE)
{
    # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
    Features = Exp = Dataset = Mean = linkColors = chrom = po = NULL


    options(stringsAsFactors = FALSE)
    figSize = 800
    segmentWidth = 25
    linePlotWidth = 90

    ##############################
    ###   networkDiagram_core.R
    ###
    ###   Authors: Michael Vacher (minor changes by Amrit :)
    ###
    ###   Parts of this code has been modified from the original OmicCircos
    ###     package obtained from:
    ###   Ying Hu Chunhua Yan <yanch@mail.nih.gov> (2015). OmicCircos:
    ###     High-quality circular visualization of omics data. v1.6.0.
    ##############################

    # check input object
    if (!is(object, "block.splsda"))
    stop("circosPlot is only available for 'block.splsda' objects")

    if (length(object$X) < 2)
    stop("This function is only available when there are more than 3 blocks
    (2 in object$X + an outcome object$Y)") # so 2 blocks in X + the outcome Y

    if (missing(cutoff))
    stop("'cutoff' is missing", call.=FALSE) # so 2 blocks in X + the outcome Y

    if(missing(color.Y))
    {
        color.Y = color.mixo(1:nlevels(object$Y))
    } else {
        if(length(color.Y) != nlevels(object$Y))
        stop("'color.Y' must be of length ", nlevels(object$Y))

    }

    if(missing(color.blocks))
    {
        color.blocks = brewer.pal(n = 12, name = 'Paired') #why 12??
    } else {
        if(length(color.blocks) != length(object$X))
        stop("'color.blocks' must be of length ", length(object$X))

        color.blocks.adj = adjustcolor(color.blocks, alpha.f = 0.5)
        #to get two shades of the same color per block

        color.blocks = c(rbind(color.blocks, color.blocks.adj))
        # to put the color next to its shaded color
    }

    if(missing(color.cor))
    {
        color.cor = c(colors()[134],  # blue, negative correlation
                colors()[128])  # pale red, positive correlation
    } else {
        if(length(color.cor) != 2)
        stop("'color.cor' must be of length 2")
    }

    X = object$X
    Y = object$Y

    #need to reorder variates and loadings to put 'Y' in last
    indY = object$indY
    object$variates = c(object$variates[-indY], object$variates[indY])
    object$loadings = c(object$loadings[-indY], object$loadings[indY])
    object$ncomp = c(object$ncomp[-indY], object$ncomp[indY])

    #check var.names
    sample.X = lapply(object$loadings[-length(object$loadings)],
    function(x){1 : nrow(x)})
    if (is.null(var.names))
    {
        var.names.list = unlist(sapply(object$loadings[
        -length(object$loadings)], rownames))
    } else if (is.list(var.names)) {
        if (length(var.names) != length(object$loadings[
        -length(object$loadings)]))
        .plotStop('var.names', sample.X)

        if(sum(sapply(1 : length(var.names), function(x){
            length(var.names[[x]]) == length(sample.X[[x]])})) !=
        length(var.names))
        .plotStop('var.names', sample.X)

        var.names.list = var.names
    } else {
        .plotStop('var.names', sample.X)
    }


    if(any(comp > min(object$ncomp)))
    {
        warning("Limitation to ",min(object$ncomp),
        " components, as determined by min(object$ncomp)")
        comp[which(comp > min(object$ncomp))] = min(object$ncomp)
    }
    comp = unique(sort(comp))


    keepA = lapply(object$loadings, function(i)
    apply(abs(i)[, comp, drop = FALSE], 1, sum) > 0)
    cord = mapply(function(x, y, keep){
        cor(x[, keep], y[, comp], use = "pairwise")
    }, x=object$X, y=object$variates[-length(object$variates)],
    keep = keepA[-length(keepA)],SIMPLIFY = FALSE)

    simMatList = vector("list", length(X))
    for(i in 1:length(cord))
    {
        for(j in 1:length(cord))
        {
            simMatList[[i]][[j]] = cord[[i]] %*% t(cord[[j]])
        }
    }
    corMat = do.call(rbind, lapply(simMatList, function(i) do.call(cbind, i)))

    ## Expression levels
    Xdat = as.data.frame(do.call(cbind, X)[, colnames(corMat)])

    AvgFeatExp0 = Xdat %>% mutate(Y = Y) %>% gather(Features, Exp, -Y) %>%
    group_by(Y, Features) %>% dplyr::summarise(Mean = mean(Exp), SD = sd(Exp))
    AvgFeatExp0$Dataset = factor(rep(names(X), unlist(lapply(cord, nrow))),
    levels = names(X))[match(AvgFeatExp0$Features,colnames(Xdat))]
    # to match Xdat that is reordered in AvgFeatExp0
    featExp = AvgFeatExp0 %>% group_by(Dataset, Y) %>% arrange(Mean)
    #Generate a circular plot (circos like) from a correlation matrix (pairwise)
    #
    # Args:
    #   corMat: the main correlation matrix.
    #         -> colnames == rownames (pairwise)  values = correlations
    #   featExp: data.frame holding the expression data.
    #   cutoff: minimum value for correlations (<threshold will be ignored)
    #   figSize: figure size
    #   segmentWidth: thickness of the segment (main circle)
    #   linePlotWidth: thickness of the line plot (showing expression data)
    #   showIntraLinks = display links intra segments



    # 1) Generate karyotype data
    chr = genChr(featExp, color.blocks = color.blocks)
    chr.names = unique(chr$chrom) # paste("chr", 1:seg.num, sep="")
    # Calculate angles and band positions
    db = segAnglePo(chr, seg=chr.names)
    db = data.frame(db)

    # 2) Generate Links
    links = genLinks(chr, corMat, threshold=cutoff)
    if (nrow(links) < 1)
    warning("Choose a lower correlation threshold to highlight
    links between datasets")

    # 3) Plot
    # Calculate parameters
    circleR = (figSize / 2.0) -  segmentWidth - linePlotWidth
    linksR = circleR - segmentWidth
    linePlotR = circleR + segmentWidth
    chrLabelsR = (figSize / 2.0)

    # replace chr$name by the ones in var.names (matching)
    # matching var.names.list with object$loadings
    ind.match = match(chr$name, unlist(sapply(object$loadings[
    -length(object$loadings)],rownames)))
    chr$name.user = unlist(var.names.list)[ind.match]

    opar1=par("mar")
    par(mar=c(2, 2, 2, 2))

    plot(c(1,figSize), c(1,figSize), type="n", axes=FALSE, xlab="",
    ylab="", main="")

    #save(list=ls(),file="temp.Rdata")
    # Plot ideogram
    drawIdeogram(R=circleR, cir=db, W=segmentWidth,  show.band.labels=TRUE,
    show.chr.labels=TRUE, chr.labels.R= chrLabelsR, chrData=chr,
    size.variables = size.variables, size.labels=size.labels,
    color.blocks = color.blocks, line = line)
    # Plot links
    if(nrow(links)>0)
    drawLinks(R=linksR, cir=db,   mapping=links,   col=linkColors,
    drawIntraChr=showIntraLinks, color.cor = color.cor)

    # Plot expression values
    cTypes = levels(Y)
    #unique(featExp[,1]) #Get the different disease/cancer types (lines)
    #lineCols = rainbow(nrow(cTypes), alpha=0.5)
    lineCols = color.Y
    #color.mixo(1:nlevels(Y))#color.mixo(match(levels(Y), levels(Y)))

    # Fixme: remove this loop and send the whole expr dframe to drawLinePlot
    if(line==TRUE)
    {
        for (i in 1:length(chr.names)){
            seg.name = gsub("chr","",chr.names[i])
            #Get data for each segment
            expr = subset(featExp,featExp$Dataset==seg.name)

            expr = dcast(expr, formula = Features ~ Y, value.var="Mean")
            ## changed PAM50 to Y
            expr = merge(expr, chr, by.x="Features", by.y="name")
            expr$po = (as.numeric(expr$chromStart) +
            as.numeric(expr$chromEnd)) / 2.0
            expr = dplyr::rename(expr, seg.name = chrom, seg.po = po)

            # Reorder columns
            cOrder = c(c(grep("seg.name", colnames(expr)),
            grep("seg.po", colnames(expr))),c(1:length(cTypes)+1))
            expr = expr[, cOrder]

            # Plot data on each sub segment
            subChr = subset(db, db$seg.name == chr.names[i] )
            drawLinePlot(R=linePlotR, cir=subChr,   W=linePlotWidth,
            lineWidth=1, mapping=expr, col=lineCols, scale=FALSE)
        }
    }
    opar=par("xpd")
    par(xpd=TRUE) # to authorise the legend to be written outside the margin,
    #       otherwise it's too small
    # Plot legend
    if(legend == TRUE)
    {
        # First legeng bottom left corner
        legend(x=5, y = (circleR/4), title="Correlations",
        c("Positive Correlation", "Negative Correlation"),
        col = color.cor, pch = 19, cex=size.legend, bty = "n")
        # Second legend bottom righ corner
        if(line==TRUE)
        legend(x=figSize-(circleR/3), y = (circleR/3), title="Expression",
        legend=levels(Y),  ## changed PAM50 to Y
        col = lineCols, pch = 19, cex=size.legend, bty = "n",ncol=ncol.legend)
        # third legend top left corner
        legend(x=figSize-(circleR/2), y = figSize, title="Correlation cut-off",
        legend=paste("r", cutoff, sep = "="),
        col = "black", cex=size.legend, bty = "n")

        legend(x=-circleR/4, y = figSize, legend=paste("Comp",
        paste(comp,collapse="-")),
        col = "black", cex=size.legend, bty = "n")
    }
    par(xpd=opar,mar=opar1)# put the previous defaut parameter for xpd
    return(invisible(corMat))
}

drawIdeogram = function(R, xc=400, yc=400, cir, W,
show.band.labels = FALSE,
show.chr.labels = FALSE, chr.labels.R = 0,
chrData,
size.variables,
size.labels,
color.blocks,
line)
{
    # Draw the main circular plot: segments, bands and labels
    chr.po    = cir
    chr.po[,1]  = gsub("chr","",chr.po[,1])
    chr.num     = nrow(chr.po)

    dat.c     = chrData
    dat.c[,1] = gsub("chr", "", dat.c[,1])

    for (chr.i in c(1:chr.num)){
        chr.s  = chr.po[chr.i,1]

        v1 = as.numeric(chr.po[chr.i,2])
        v2 = as.numeric(chr.po[chr.i,3])
        v3 = as.numeric(chr.po[chr.i,6])
        v4 = as.numeric(chr.po[chr.i,7])

        dat.v = subset(dat.c, dat.c[,1]==chr.s)
        dat.v = dat.v[order(as.numeric(dat.v[,2])),]
        for (i in 1:nrow(dat.v)){

            #col.v = which(colors()==dat.v[i,5])  #get color index
            #col = colors()[col.v]
            dark.clear = color.blocks#brewer.pal(n = 12, name = 'Paired')
            col.v = which(dark.clear==dat.v[i,5])  #get color index
            col = dark.clear[col.v]

            w1 = scale.v(as.numeric(dat.v[i,2]), v1, v2, v3, v4)
            w2 = scale.v(as.numeric(dat.v[i,3]), v1, v2, v3, v4)

            draw.arc.s(xc, yc, R, w1, w2, col=col, lwd=W)

            if (show.band.labels){
                band.text = as.character(dat.v[i,"name.user"])

                band.po = ((w1+w2)/2)# - ((w2-w1)/3) #position around the circle
                # print(c(band.po, w1, w2, (w2-w1)/3))
                band.po.in = R-(W/3.0) #position on the band (middle)
                draw.text.rt(xc, yc,band.po.in  , band.po , band.text ,
                cex = size.variables, segmentWidth = W, side="in" )
            }
        } #End for row
        if (show.chr.labels){
            w.m = (v1+v2)/2
            chr.t = gsub("chr", "", chr.s)
            if(line == TRUE)
            {
                draw.text.rt(xc, yc, chr.labels.R, w.m, chr.t, cex=size.labels,
                segmentWidth = W, parallel=TRUE)
            } else {
                #put the labels closer to the circle
                draw.text.rt(xc, xc, chr.labels.R, w.m, chr.t, cex=size.labels,
                segmentWidth = 75, parallel=TRUE)
            }
        }
    } #End for
}

drawLinks = function(R, xc=400, yc=400, cir, W,
mapping=mapping,
lineWidth=1, col=rainbow(10, alpha=0.8)[7],  drawIntraChr=FALSE,
color.cor = color.cor)
{
    # Draw the links (computed correlation) between features
    chr.po    = cir
    chr.po[,1]  = gsub("chr","",chr.po[,1])
    chr.num     = nrow(chr.po)

    chr.po[,4] = gsub("chr", "", chr.po[,4])
    dat.in = mapping
    dat.in[,1] = gsub("chr", "", dat.in[,1])
    dat.in[,4] = gsub("chr", "", dat.in[,4])

    dat    = dat.in

    for (i in 1:nrow(dat)){
        chr1.s   = dat[i,1]
        chr2.s   = dat[i,4]
        po1      = dat[i,2]
        po2      = dat[i,5]

        chr1     = which(chr.po[,1]==chr1.s)
        chr2     = which(chr.po[,1]==chr2.s)

        v1 = as.numeric(chr.po[chr1,2])
        v2 = as.numeric(chr.po[chr1,3])
        v3 = as.numeric(chr.po[chr1,6])
        v4 = as.numeric(chr.po[chr1,7])

        w1 = scale.v(as.numeric(po1), v1, v2, v3, v4)

        v1 = as.numeric(chr.po[chr2,2])
        v2 = as.numeric(chr.po[chr2,3])
        v3 = as.numeric(chr.po[chr2,6])
        v4 = as.numeric(chr.po[chr2,7])

        w2 = scale.v(as.numeric(po2), v1, v2, v3, v4)
        # Set the link width depending on the correlation coefficient
        lwd = abs(as.numeric(dat[i,7]))

        # Set link color
        if (as.numeric(dat[i,7]) < 0.0){
            linkCol = color.cor[2]#colors()[128]  #pale red
        } else {
            linkCol = color.cor[1]#colors()[134]  # blue
        }
        linkCol = add.alpha(linkCol, alpha=0.4)

        if (chr1 == chr2){
            if (drawIntraChr == TRUE){
                draw.link(xc, yc, R, w1, w2, col=linkCol, lwd=lineWidth)
            }
        } else {
            draw.link(xc, yc, R, w1, w2, col=linkCol, lwd=lineWidth)
        }
    } ### End for
}

drawLinePlot = function(mapping=mapping, xc=400, yc=400, col.v=3,
R, cir,   W, col='black', scale=FALSE, lineWidth=1,
background.lines=FALSE,axis.width=1)
{
    # Generate a linear plot around the main ideogram.
    #
    # fixme: the function writes the same line multiple times
    # it needs to be called only once and process the data/segment
    # separately
    chr.po    = cir
    chr.po[,1]  = gsub("chr","",chr.po[,1])
    chr.num     = nrow(chr.po)

    dat.in   = mapping
    dat.in[,1] = gsub("chr", "", dat.in[,1])

    # data set for the chromosome
    for (chr.i in 1:chr.num){
        chr.s = chr.po[chr.i,1]
        chr.s = gsub("chr","",chr.s)
        dat   = subset(dat.in, dat.in[,1]==chr.s)
        dat   = dat[order(as.numeric(dat[,2])),]
        v1 = as.numeric(chr.po[chr.i,2])
        v2 = as.numeric(chr.po[chr.i,3])
        v3 = as.numeric(chr.po[chr.i,6])
        v4 = as.numeric(chr.po[chr.i,7])

        # background line
        if (background.lines){
            draw.arc.pg(xc, yc, v1, v2, R, R+W-5, col=colors()[245])
        } else {
            draw.arc.s(xc, yc, R, v1, v2, col=colors()[245], lwd=axis.width)
        }
    }

    my.R1 = R + W/5
    my.R2 = R + W - W/5

    ## for the matrix colors
    num.col = ncol(dat[,col.v:ncol(dat)])
    num.row = nrow(dat.in)

    if (length(col) == num.col){
        colors = col
    } else {
        colors = rainbow(num.col, alpha=0.5)
    }

    for (chr.i in 1:chr.num){
        chr.s = chr.po[chr.i,1]
        chr.s = gsub("chr","",chr.s)

        dat   = subset(dat.in, dat.in[,1]==chr.s)
        dat   = dat[order(as.numeric(dat[,2])),]
        #print(head(dat))
        dat.i   = c(col.v:ncol(dat))
        dat.m   = dat.in[,dat.i]
        dat.m   = as.matrix(dat.m)
        dat.min = min(as.numeric(dat.m), na.rm=TRUE)
        dat.max = max(as.numeric(dat.m), na.rm=TRUE)

        v1 = as.numeric(chr.po[chr.i,2])
        v2 = as.numeric(chr.po[chr.i,3])
        v3 = as.numeric(chr.po[chr.i,6])
        v4 = as.numeric(chr.po[chr.i,7])

        col.i = 0

        for (j in col.v:ncol(dat)){
            col.i = col.i + 1
            col   = colors[col.i]

            my.v      = as.numeric(dat[1,j])
            dat.i.old = my.v
            v.old   = scale.v(my.v, my.R1, my.R2, dat.min, dat.max)

            po      = as.numeric(dat[1,2])
            w.from  = scale.v(po, v1, v2, v3, v4)


            for (k in 1:nrow(dat)){

                dat.i = as.numeric(dat[k,j])

                if (is.na(dat.i)){
                    next
                }

                v    = scale.v(dat.i, my.R1, my.R2, dat.min, dat.max)
                w.to = scale.v(as.numeric(dat[k,2]), v1, v2, v3, v4)

                if (w.from > 0){
                    draw.line3(xc, yc, w.from, w.to, v.old, v, col=col,
                    lwd=lineWidth)
                }

                dat.i.old = dat.i
                w.from    = w.to
                v.old     = v
            } # end the row
        }   # end the col
        if (scale){
            do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5))
        }
    }     # end the chr/segment
}

genChr =function (expr, bandWidth = 1.0, color.blocks)
{
    # Generate the segments and calculate the
    # unique positions of the bands
    #
    # Args:
    #   expr : dataframe containing the features expression
    # example: colnames(concatFeatExp) "PAM50"    "Features" "Mean"     "SD"
    # "Dataset"
    #   bandWidth: thickness of each band
    #
    # Return:
    #   a data.frame that can be used with segAnglePo

    # expr can contains expression data for multiple diseases
    # here, we only use the Chrom and Dataset column and remove duplicates
    keeps = c("Features","Dataset")
    expr = expr[keeps]
    expr = unique(expr)
    chrLengths = data.frame(table(expr$Dataset))
    rownames(chrLengths) = chrLengths[,1]
    chrLengths[,1] = NULL
    colnames(chrLengths) = c( "Freq")
    chrLengths[, "Count"] = rep(0, nrow(chrLengths))

    # Last column contains the bands' color
    # Create a color scheme
    #dark = c("brown3","darkgoldenrod","antiquewhite3","steelblue3")
    #clear = c("brown1","darkgoldenrod1","antiquewhite1","steelblue1")
    dark.clear = color.blocks#brewer.pal(n = 12, name = 'Paired')
    dark = dark.clear[seq(2, 12, by = 2)]
    clear = dark.clear[seq(1, 12, by = 2)]
    chrColScheme = data.frame(dark, clear)
    n_datasets = length(unique(expr$Dataset))
    chrColScheme = chrColScheme[c(1:n_datasets),]
    rownames(chrColScheme) = levels(factor(expr$Dataset))#alphabetical order

    seg.out = c()
    for (i in 1:nrow(expr)){
        chrName    = paste("chr", as.character(expr[i,'Dataset'][[1]]), sep="")
        dType = as.character(expr[i,'Dataset'][[1]])
        pStart = chrLengths[dType,'Count'] * bandWidth
        pStop = chrLengths[dType,'Count'] * bandWidth + bandWidth
        chrLengths[dType,'Count'] = chrLengths[dType,'Count'] + 1
        fName = as.character(as.matrix(expr[i,'Features']))
        # added as.character() by amrit
        # Assign colors
        if (chrLengths[dType,'Count'] %% 2 == 0){
            chrCol = chrColScheme[dType,]$clear
        } else{
            chrCol = chrColScheme[dType,]$dark
        }
        seg.out = rbind(seg.out, c(chrName, pStart, pStop,  fName, chrCol))
    }

    # Use the same names than in omicCircos
    colnames(seg.out) = c("chrom", "chromStart", "chromEnd", "name", "color")
    seg.out = as.data.frame(seg.out)

    return(seg.out)
}

genLinks = function(chr, corMat, threshold)
{

    # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
    Var1=Var2=chrom=NULL


    # Generates the links corresponding to pairwise correlations
    #
    # Args:
    #   chr: ideogram structure (generated from genChr)
    #   corMat: main correlation matrix
    #   threshold: minimum correlation value
    #
    # Return:
    #   the links data (see omicsCircos doc)
    #
    linkList = c()
    # Remove matrix diagonal and the the mat in a list
    linkList = subset(melt(corMat), Var1!=Var2 )
    # Remove links below the threshold
    linkList = subset(linkList, abs(linkList$value) >= threshold)

    #First merge
    linkList = dplyr::rename(linkList, feat1=Var1, feat2=Var2)
    # CHANGED BY AMRIT
    linkList = merge(linkList, chr, by.x="feat1", by.y="name")
    # Set the position in the middle of the band
    linkList$po1 = (as.numeric(linkList$chromStart)
    + as.numeric(linkList$chromEnd)) / 2.0
    linkList = dplyr::rename(linkList, chr1=chrom)   # CHANGED BY AMRIT
    keeps = c("feat1","feat2","value","chr1","po1")
    linkList = linkList[keeps]

    #Second merge
    linkList = merge(linkList, chr, by.x="feat2", by.y="name")
    linkList$po2 = (as.numeric(linkList$chromStart)
    + as.numeric(linkList$chromEnd)) / 2.0
    linkList = dplyr::rename(linkList, chr2=chrom)   # CHANGED BY AMRIT
    keeps = c("chr1","po1","feat1","chr2","po2","feat2","value")
    linkList = linkList[keeps]

    return(linkList)
}

bezierCurve = function(x, y, n=10)  {
    outx = NULL
    outy = NULL
    i = 1
    for (t in seq(0, 1, length.out=n))		{
        b = bez(x, y, t)
        outx[i] = b$x
        outy[i] = b$y
        i = i+1
    }
    return (list(x=outx, y=outy))
}

##
bez = function(x, y, t)	{
    outx = 0
    outy = 0
    n = length(x)-1
    for (i in 0:n)		{
        outx = outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
        outy = outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
    }
    return (list(x=outx, y=outy))
}

###########################################
# one value : from a to b
scale.v = function(v, a, b, min.v, max.v) {
    v = v-min.v
    v = v/(max.v-min.v)
    v = v*(b-a)
    v+a
}

### draw.link
draw.link = function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
    # for translocation
    w3  = (w1+w2)/2
    w1  = w1/360*2*pi
    w2  = w2/360*2*pi
    w3  = w3/360*2*pi
    x0  = xc+r*cos(w1)
    y0  = yc-r*sin(w1)
    x1  = xc+r*cos(w2)
    y1  = yc-r*sin(w2)
    x = c(x0,xc,xc,x1)
    y = c(y0,yc,yc,y1)
    points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}

### draw.link2
draw.link2 = function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
    # for translocation
    w3  = (w1+w2)/2
    w1  = w1/360*2*pi
    w2  = w2/360*2*pi
    w3  = w3/360*2*pi
    x2  = xc+r/2*cos(w3)
    y2  = yc-r/2*sin(w3)
    x0  = xc+r*cos(w1)
    y0  = yc-r*sin(w1)
    x1  = xc+r*cos(w2)
    y1  = yc-r*sin(w2)
    x = c(x0, x2, x2, x1)
    y = c(y0, y2, y2, y1)
    points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}
###

###
draw.link.pg = function(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd) {
    w1 = w1.1
    w2 = w2.2
    w3  = (w1+w2)/2
    w1  = w1/360*2*pi
    w2  = w2/360*2*pi
    w3  = w3/360*2*pi
    x0  = xc+r*cos(w1)
    y0  = yc-r*sin(w1)
    x1  = xc+r*cos(w2)
    y1  = yc-r*sin(w2)
    x = c(x0,xc,xc,x1)
    y = c(y0,yc,yc,y1)
    bc1 = bezierCurve(x,y,60)

    ang.d = abs(w1.1-w1.2)
    pix.n = ang.d * 10
    if (pix.n < 10){
        pix.n = 10
    }

    ang.seq = rev(seq(w1.1,w1.2,length.out=pix.n))
    ang.seq = ang.seq/360*2*pi

    fan.1.x = xc + cos(ang.seq) * r
    fan.1.y = yc - sin(ang.seq) * r

    ######################################################
    w1 = w1.2
    w2 = w2.1
    w3  = (w1+w2)/2
    w1  = w1/360*2*pi
    w2  = w2/360*2*pi
    w3  = w3/360*2*pi
    x0  = xc+r*cos(w1)
    y0  = yc-r*sin(w1)
    x1  = xc+r*cos(w2)
    y1  = yc-r*sin(w2)
    x = c(x0,xc,xc,x1)
    y = c(y0,yc,yc,y1)
    bc2 = bezierCurve(x,y,60)

    ang.d = abs(w2.1-w2.2)
    pix.n = ang.d * 10
    if (pix.n < 10){
        pix.n = 10
    }

    ang.seq = rev(seq(w2.1,w2.2,length.out=pix.n))
    ang.seq = ang.seq/360*2*pi

    fan.2.x = xc + cos(ang.seq) * r
    fan.2.y = yc - sin(ang.seq) * r

    polygon(c(bc1$x, fan.2.x, rev(bc2$x), rev(fan.1.x)),
    c(bc1$y, fan.2.y, rev(bc2$y), rev(fan.1.y)),
    fillOddEven=TRUE, border=col, col=col, lwd=lwd)
}

###
draw.point.w = function(xc, yc, r, w, col=col, cex=cex){
    w = w/360*2*pi
    x = xc+r*cos(w)
    y = yc-r*sin(w)
    points(x, y, pch=20, col=col, cex=cex)
}

###
draw.text.w = function(xc, yc, r, w, n, col="black", cex=1){
    w = w%%360
    w = w/360*2*pi
    x = xc+r*cos(w)
    y = yc-r*sin(w)
    text(x,y,labels=n, col=col, cex=cex)
}

###
draw.text.rt = function(xc, yc, r, w, n, col="black", cex=1, side="out",
segmentWidth=20, parallel=FALSE){
    w     = w%%360
    the.o = w

    the.w = 360-w
    w     = w/360*2*pi
    x     = xc+r*cos(w)
    y     = yc-r*sin(w)


    num2  = (segmentWidth*2)/2.0
    b = the.w
    if (side=="out"){
        if (the.w <= 90 ){
            the.pos = 4
            if (parallel == TRUE){
                the.w = the.w -90 # 180
            }
        } else if (the.w > 90 & the.w <= 180) {
            if (parallel == TRUE){
                the.w = the.w -90 # 180
            }
            else {
                the.w = the.w + 180
            }
            the.pos = 2
        } else if (the.w > 180 & the.w <= 270){
            the.w = the.w%%180
            if (parallel == TRUE){
                the.w = the.w -90 # 180
            }
            the.pos = 2
        } else if (the.w > 270 & the.w <= 360){
            the.pos = 4
            if (parallel == TRUE){
                the.w = the.w + 90
            }
        }

        if (the.pos==2){
            x = x+num2
        }
        if (the.pos==4){
            x = x-num2
        }
    }

    if (side=="in"){
        if (the.w <= 90 ){
            the.pos = 4
        } else if (the.w > 90 & the.w <= 180) {
            the.w = the.w + 180
            the.pos = 2
        } else if (the.w > 180 & the.w <= 270){
            the.w = the.w%%180
            the.pos = 2
        } else if (the.w > 270 & the.w <= 360){
            the.pos = 4
        }

        if (the.pos==2){
            x = x+segmentWidth
        }
        if (the.pos==4){
            x = x-segmentWidth
        }
    }


    text(x, y, adj=0, offset=1, labels=n, srt=the.w,
    pos=the.pos, col=col, cex=cex)
}

###strokeLine2
draw.line = function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
    w  = (w/360)*2*pi
    x1 = xc+l1*cos(w)
    y1 = yc-l1*sin(w)
    x2 = xc+l2*cos(w)
    y2 = yc-l2*sin(w)
    segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend)
}

###strokeLine3
draw.line2 = function (xc, yc, w, r, l, col=col, lwd=lwd){
    line_w   = l
    theangle = w
    l1       = r
    theangle = (theangle/360)*2*pi
    x0       = xc+l1*cos(theangle)
    y0       = yc+l1*sin(theangle)
    w1       = 45/360*2*pi
    x1 = xc + sin(w1) * (x0)
    y1 = yc + cos(w1) * (y0)
    x2 = xc - sin(w1) * (x0)
    y2 = yc - cos(w1) * (y0)
    segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt")
}

###strokeLine by two angles
draw.line3 = function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd){
    theangle1 = w1
    theangle2 = w2
    l1        = r1
    l2        = r2

    theangle1 = (theangle1/360)*2*pi
    x1        = xc+l1*cos(theangle1)
    y1        = yc-l1*sin(theangle1)

    theangle2 = (theangle2/360)*2*pi
    x2        = xc+l2*cos(theangle2)
    y2        = yc-l2*sin(theangle2)

    segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt")
}

### plot fan or sector that likes a piece of doughnut (plotFan)
draw.arc.pg = function (xc, yc,
w1, w2, r1, r2, col="lightblue", border="lightblue", lwd=0.01
){

    ang.d = abs(w1-w2)
    pix.n = ang.d * 10
    if (pix.n < 10){
        pix.n = 10
    }

    ang.seq = rev(seq(w1,w2,length.out=pix.n))
    ang.seq = ang.seq/360*2*pi

    fan.i.x = xc + cos(ang.seq) * r1
    fan.i.y = yc - sin(ang.seq) * r1


    fan.o.x = xc + cos(ang.seq) * r2
    fan.o.y = yc - sin(ang.seq) * r2

    polygon(c(rev(fan.i.x), fan.o.x ), c(rev(fan.i.y), fan.o.y),
    fillOddEven=TRUE, border=border, col=col, lwd=lwd, lend=1)

}

draw.arc.s = function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
    # Draw circular arcs for the main ideogram)
    # s = simple
    # r = radius
    ang.d = abs(w1-w2)
    pix.n = ang.d * 5
    if (pix.n < 2){
        pix.n = 2
    }

    ang.seq = rev(seq(w1,w2,length.out=pix.n))
    ang.seq = ang.seq/360*2*pi

    fan.i.x = xc + cos(ang.seq) * r
    fan.i.y = yc - sin(ang.seq) * r
    ## lend=0(round)  #lend=1(butt)  lend=2(square)
    lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend)
    #points(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend)
}

#########################################
## segment to angle and position
## segAnglePo
#########################################

# get angle if given seg and position
# seg should be ordered by user
segAnglePo = function (seg.dat=seg.dat, seg=seg, angle.start=angle.start,
angle.end=angle.end){

    if (missing(angle.start)){
        angle.start = 0
    }
    if (missing(angle.end)){
        angle.end = 360
    }
    ## check data.frame?
    colnames(seg.dat) = c("seg.name","seg.Start","seg.End","name","gieStain")

    ## get length of the segomosomes
    seg.l   = c()
    seg.min = c()
    seg.sum = c()
    seg.s   = 0
    seg.num   = length(seg)
    seg.names = seg

    ########################################################
    ########################################################
    for (i in 1:seg.num){
        seg.n = seg.names[[i]]

        dat.m = subset(seg.dat, seg.dat[,1]==seg.n)
        seg.full.l = max(as.numeric(dat.m[,"seg.End"]))
        seg.full.m = min(as.numeric(dat.m[,"seg.Start"]))
        seg.l      = cbind(seg.l, seg.full.l)
        seg.min    = cbind(seg.min, seg.full.m)
        seg.s      = seg.s + seg.full.l
        seg.sum    = cbind(seg.sum, seg.s)
    }

    ## initial parameters
    gap.angle.size = 2
    seg.angle.from = angle.start + 270

    seg.full.l  = sum(as.numeric(seg.l))
    angle.range = angle.end - angle.start
    cir.angle.r = (angle.range - seg.num * gap.angle.size)/seg.full.l

    out.s     = c()
    l.old     = 0
    gap.angle = 0
    for (i in 1:seg.num){
        seg.n = seg.names[[i]]
        dat.m = subset(seg.dat, seg.dat[,1]==seg.n)
        len   = seg.sum[i]
        w1    = cir.angle.r*l.old + gap.angle
        w2    = cir.angle.r*len   + gap.angle
        out.s     = rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from,
        l.old, len, seg.min[i], seg.l[i]))
        gap.angle = gap.angle + gap.angle.size
        l.old     = len
    }

    colnames(out.s) = c("seg.name","angle.start", "angle.end", "seg.sum.start",
    "seg.sum.end","seg.start", "seg.end")
    return(out.s)
}


### start do.scale
do.scale = function(xc=xc, yc=yc, dat.min=dat.min, dat.max=dat.max,
R=R, W=W, s.n=1, col="blue"){
    dat.m   = round((dat.min+dat.max)/2, s.n)
    dat.min = round(dat.min, s.n)
    dat.max = round(dat.max, s.n)
    y1      = yc + R
    y2      = yc + R + W/2
    y3      = yc + R + W
    x1      = xc - W/20
    x2      = x1 - (W/20)*1.2
    x3      = x1 - (W/20)*3
    segments(x1, y1, x1, y3, lwd=0.01, col=col)
    segments(x1, y1, x2, y1, lwd=0.01, col=col)
    segments(x1, y2, x2, y2, lwd=0.01, col=col)
    segments(x1, y3, x2, y3, lwd=0.01, col=col)
    text(x3, y1, dat.min, cex=0.2, col=col)
    text(x3, y2, dat.m,   cex=0.2, col=col)
    text(x3, y3, dat.max, cex=0.2, col=col)
}
### end do.scale

## Add an alpha value to a colour
add.alpha = function(col, alpha=1){
    if(missing(col))
    stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
    function(x)
    rgb(x[1], x[2], x[3], alpha=alpha))
}
###-------------------------------------------------------------------------


