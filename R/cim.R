################################################################################
# Authors:
#   Ignacio Gonzalez,
#   Francois Bartolo,
#   Kim-Anh Le Cao,
# created: 2013
# last modified: 2015
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
################################################################################


# -------------------------------------------- #
# CIM for objects "pca","spca","ipca","sipca","mlsplsda","splsda","plsda","rcc",
#   "pls","spls","mlspls"
# -------------------------------------------- #
cim <-
    function(mat=NULL,
             color = NULL,
             row.names = TRUE,
             col.names = TRUE,
             row.sideColors = NULL,
             col.sideColors = NULL,
             row.cex = NULL,
             col.cex = NULL,
             threshold = 0,
             cluster = "both",
             dist.method = c("euclidean", "euclidean"),
             clust.method = c("complete", "complete"),
             cut.tree = c(0, 0),
             transpose = FALSE,
             symkey = TRUE,
             keysize = c(1, 1),
             keysize.label = 1,
             zoom = FALSE,
             title = NULL,
             xlab = NULL,
             ylab = NULL,
             margins = c(5, 5),
             lhei = NULL,
             lwid = NULL,
             comp = NULL,
             center = TRUE,
             scale = FALSE,
             mapping = "XY",
             legend = NULL,
             save = NULL,
             name.save = NULL)

    {
        
    class.object <- class(mat)
    
    
    
    #-- checking general input parameters -------------------------------------
    #--------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
    error = function(e) e)
    
    if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)

    #-- color
    if (is.null(color))
    color = color.spectral(25)
    
    
    #-- cluster
    choices = c("both", "row", "column", "none")
    cluster = choices[pmatch(cluster, choices)]
    
    if (is.na(cluster))
    stop("'cluster' should be one of 'both', 'row', 'column' or 'none'.",
    call. = FALSE)
    
    #-- cluster method
    if (!is.character(clust.method) | length(as.vector(clust.method)) != 2)
    stop("'clust.method' must be a character vector of length 2.",
    call. = FALSE)
    
    choices = c("ward.D", "single", "complete", "average", "mcquitty",
    "median", "centroid")
    clust.method = choices[c(pmatch(clust.method[1], choices),
    pmatch(clust.method[2], choices))]
    
    if (any(is.na(clust.method)))
    stop("invalid clustering method.", call. = FALSE)
    
    #-- distance method
    if (!is.character(dist.method) | length(as.vector(dist.method)) != 2)
    stop("'dist.method' must be a character vector of length 2.", call. = FALSE)
    
    choices = c("euclidean", "correlation", "maximum", "manhattan",
    "canberra", "binary", "minkowski")
    dist.method = choices[c(pmatch(dist.method[1], choices),
    pmatch(dist.method[2], choices))]
    
    if (any(is.na(dist.method)))
    stop("invalid distance method.", call. = FALSE)
    
    
    #-- checking general input arguments -------------------------------------
    #-------------------------------------------------------------------------#
    
    #-- color
    if (any(!sapply(color, function(color) {
        tryCatch(is.matrix(col2rgb(color)), error = function(e) FALSE) })))
    stop("'color' must be a character vector of recognized colors.",
    call. = FALSE)
    
    
    #-- row.sideColors
    if (any(!sapply(row.sideColors, function(row.sideColors) {
        tryCatch(is.matrix(col2rgb(row.sideColors)),
        error = function(e) FALSE) })))
    stop("color names for vertical side bar must be a character vector
    of recognized colors.",
    call. = FALSE)
    
    
    #-- col.sideColors
    if (any(!sapply(col.sideColors, function(col.sideColors) {
        tryCatch(is.matrix(col2rgb(col.sideColors)),
        error = function(e) FALSE) })))
    stop("color names for horizontal side bar must be a character vector
    of recognized colors.",
    call. = FALSE)
    
    
    #-- row.cex
    if (!is.null(row.cex)) {
        if (!is.numeric(row.cex) || length(row.cex) != 1)
        stop("'row.cex' must be a numerical value.", call. = FALSE)
    }
    
    #-- col.cex
    if (!is.null(col.cex)) {
        if (!is.numeric(col.cex) || length(col.cex) != 1)
        stop("'col.cex' must be a numerical value.", call. = FALSE)
    }
    
    #-- transpose
    if (!is.logical(transpose))
    stop("'transpose' must be a logical constant (TRUE or FALSE).",
    call. = FALSE)
    
    #-- cut.tree
    if (!is.numeric(cut.tree) || length(cut.tree) != 2)
    stop("'cut.tree' must be a numeric vector of length 2.",
    call. = FALSE)
    else {
        if (!(all(0 <= cut.tree & cut.tree <= 1)))
        stop("Components of 'cut.tree' must be between 0 and 1.",
        call. = FALSE)
    }
    
    #-- keysize
    if (length(keysize) != 2 || any(!is.finite(keysize)))
    stop("'keysize' must be a numeric vector of length 2.",
    call. = FALSE)
    
    #-- keysize.label
    if (length(keysize.label) != 1 || any(!is.finite(keysize)))
    stop("'keysize' must be a numeric vector of length 1.",
    call. = FALSE)
    
    #-- zoom
    if (!is.logical(zoom))
    stop("'zoom' must be a logical constant (TRUE or FALSE).",
    call. = FALSE)
    
    #-- margins
    if (!is.numeric(margins) || length(margins) != 2)
    stop("'margins' must be a numeric vector of length 2.",
    call. = FALSE)
    
    #-- symkey
    if (!is.logical(symkey))
    stop("'symkey' must be a logical constant (TRUE or FALSE).",
    call. = FALSE)
    
    #-- lhei
    if (!is.null(lhei)) {
        if (is.null(col.sideColors)) {
            if (length(lhei) != 2 | !is.numeric(lhei) | any(is.na(lhei)))
            stop("'lhei' must be a numeric vector of length 2.",
            call. = FALSE)
        }
        else {
            if (length(lhei) != 3 | !is.numeric(lhei) | any(is.na(lhei)))
            stop("'lhei' must be a numeric vector of length 3.",
            call. = FALSE)
        }
    }
    
    #-- lwid
    if (!is.null(lwid)) {
        if (is.null(row.sideColors)) {
            if (length(lwid) != 2 | !is.numeric(lwid) | any(is.na(lwid)))
            stop("'lwid' must be a numeric vector of length 2.",
            call. = FALSE)
        }
        else {
            if (length(lwid) != 3 | !is.numeric(lwid) | any(is.na(lwid)))
            stop("'lwid' must be a numeric vector of length 3.",
            call. = FALSE)
        }
    }
    
    #-- xlab
    xlab = as.graphicsAnnot(xlab)
    
    #-- ylab
    ylab = as.graphicsAnnot(ylab)
    
    #-- title
    title = as.graphicsAnnot(title)
    
    #-- threshold correlation
    if (!is.numeric(threshold) | (threshold > 1) | (threshold < 0))
    stop("The value taken by 'threshold' must be between 0 and 1",
    call. = FALSE)
    
    #-- save
    if (!is.null(save)){
        if (! save %in% c("jpeg","tiff","png","pdf"))
        stop("'save' must be one of 'jpeg', 'png', 'tiff' or 'pdf'.",
        call. = FALSE)
    }
    
    #-- name.save
    if (!is.null(name.save)){
        if (! is.character(name.save) || length(name.save) > 1)
        stop("'name.save' must be a character.", call. = FALSE)
    } else {
        if (!is.null(save))
        name.save = paste0("cim_",gsub(".", "_", deparse(substitute(mat)),
        fixed = TRUE))
    }
    
    #-- end checking --#
    #------------------#
    if (!is.null(save)){
        
        while (dev.cur()>2)
        dev.off()
        
        if (save == "jpeg")
        jpeg(filename = paste0(name.save,".jpeg"), res = 600, width = 4000,
        height = 4000)
        if (save == "png")
        jpeg(filename = paste0(name.save,".png"), res = 600, width = 4000,
        height = 4000)
        if (save == "tiff")
        tiff(filename = paste0(name.save,".tiff"), res = 600, width = 4000,
        height = 4000)
        if (save == "pdf")
        pdf(file = paste0(name.save,".pdf"))
        
    }
    
    
    object.pca=c("pca","spca","ipca","sipca","mixo_mlsplsda","mixo_splsda",
    "mixo_plsda")
    object.rcc=c("rcc")
    object.pls=c("mixo_pls","mixo_spls","mixo_mlspls")
    object.list=c("pca","spca","ipca","sipca","mixo_mlsplsda","mixo_splsda",
    "mixo_plsda", "rcc","mixo_pls","mixo_spls","mixo_mlspls")
    
    if (any(class.object == "block.splsda"))
    stop("Please call the 'cimDiablo' function on your 'block.splsda' object",
    call. = FALSE)
    
    
    if (!any(class.object %in% c(object.list,"matrix")))
    stop("'mat' has to be a matrix or one of the following object: ",
    paste(object.list, collapse =", "), ".", call. = FALSE)
    
    #-- if mixOmics class
    if(any(class.object  %in%  object.list))
    {
        #-- general checks  -------------
        p = ncol(mat$X)
        q = ncol(mat$Y)
        n = nrow(mat$X)
        ncomp = mat$ncomp
        #-- comp
        if (is.null(comp)) {
            comp = 1:mat$ncomp
        }
        if (length(comp) > 1) {
            comp = unique(comp)
            if (!is.numeric(comp) || any(comp < 1))
                stop("invalid vector for 'comp'.", call. = FALSE)
            if (any(comp > ncomp))
                stop("the elements of 'comp' must be smaller or equal than ",
                     ncomp, ".", call. = FALSE)
        }
        
        if (length(comp) == 1) {
            if (is.null(comp) || !is.numeric(comp) || comp <= 0 || comp > ncomp)
                stop("invalid value for 'comp'.", call. = FALSE)
            comp <- c(comp, comp)
        }
        
        comp <- round(comp)
        
        # if object is a pls or spls with a univariate Y, or multivariate but
        # only 1 variable selected on all comp, we only plot a heatmap of X
        if(any(class.object %in%  object.pls))
        {
            temp = apply(mat$loadings$Y, 2,
                         function(x) {
                             which(x != 0, arr.ind = TRUE)
                         })
            # gives which variables are selected
            num.variable.selected.Y = table(unlist(temp))
            
            if (length(num.variable.selected.Y) == 1)
                #only one variable in Y to plot, will raise trouble so we
                # switch from (s)pls to (s)plsda
                class.object = "mixo_splsda"
        }
        #-- if c("mixo_pls","mixo_spls","mixo_mlspls") or pls with univarite Y ----
        ## or multivariate but only one Y kept in sparse model
        
        if (!any(class.object  %in%  object.pca)) {
            #-- mapping
            choices = c("XY", "X", "Y")
            mapping = choices[pmatch(mapping, choices)]
            
            if (is.na(mapping))
                stop("'mapping' should be one of 'XY', 'X' or 'Y'.", call. = FALSE)
            
            if (mapping == "XY")
            {
                if (is.logical(row.names))
                {
                    if (isTRUE(row.names))
                    row.names = mat$names$colnames$X
                    else
                    row.names = rep("", p)
                } else {
                    if (length(row.names) != p)
                    stop("'row.names' must be a character vector of length ",
                    p, ".", call. = FALSE)
                }
                
                if (is.logical(col.names))
                {
                    if (isTRUE(col.names))
                    col.names = mat$names$colnames$Y
                    else
                    col.names = rep("", q)
                } else {
                    if (length(col.names) != q)
                    stop("'col.names' must be a character vector of length ",
                    q, ".", call. = FALSE)
                }
                
                if (!is.null(row.sideColors))
                {
                    row.sideColors = as.matrix(row.sideColors)
                    if (nrow(row.sideColors) != p)
                    stop("'row.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ", p, ".",
                    call. = FALSE)
                }
                if (!is.null(col.sideColors))
                {
                    col.sideColors = as.matrix(col.sideColors)
                    if (nrow(col.sideColors) != q)
                    stop("'col.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ", q, ".",
                    call. = FALSE)
                }
            }
            
            
            if (mapping == "X")
            {
                if (is.logical(row.names))
                {
                    if (isTRUE(row.names)) {
                        row.names = mat$names$sample
                    }
                    else
                    row.names = rep("", n)
                } else {
                    if (length(row.names) != n)
                    stop("'row.names' must be a character vector of length ",
                    n, ".",
                    call. = FALSE)
                }
                if (is.logical(col.names))
                {
                    if (isTRUE(col.names))
                    col.names = mat$names$colnames$X
                    else
                    col.names = rep("", p)
                } else {
                    if (length(col.names) != p)
                    stop("'col.names' must be a character vector of length ",
                    p, ".",
                    call. = FALSE)
                }
                
                if (!is.null(row.sideColors))
                {
                    row.sideColors = as.matrix(row.sideColors)
                    if (nrow(row.sideColors) != n)
                    stop("'row.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ", n, ".",
                    call. = FALSE)
                }
                if (!is.null(col.sideColors))
                {
                    col.sideColors = as.matrix(col.sideColors)
                    if (nrow(col.sideColors) != p)
                    stop("'col.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ", p, ".",
                    call. = FALSE)
                }
            }
            
            if (mapping == "Y")
            {
                if (is.logical(row.names))
                {
                    if (isTRUE(row.names))
                    {
                        if (any(class.object %in% object.rcc))
                        row.names = mat$names$sample
                        else row.names = mat$names$sample
                    }
                    else
                    row.names = rep("", n)
                } else {
                    if (length(row.names) != n)
                    stop("'row.names' must be a character vector of length ",
                    n, ".",
                    call. = FALSE)
                }
                if (is.logical(col.names))
                {
                    if (isTRUE(col.names))
                    col.names = mat$names$colnames$Y
                    else
                    col.names = rep("", q)
                } else {
                    if (length(col.names) != q)
                    stop("'col.names' must be a character vector of length ",
                    q, ".",
                    call. = FALSE)
                }
                
                if (!is.null(row.sideColors))
                {
                    row.sideColors = as.matrix(row.sideColors)
                    if (nrow(row.sideColors) != n)
                    stop("'row.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ", n, ".",
                    call. = FALSE)
                }
                if (!is.null(col.sideColors))
                {
                    col.sideColors = as.matrix(col.sideColors)
                    if (nrow(col.sideColors) != q)
                    stop("'col.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ", q, ".",
                    call. = FALSE)
                }
            }
        }
        
        #-- if NOT c("mixo_pls","mixo_spls","mixo_mlspls") or pls with univarite Y ----
        if(any(class.object %in%  object.pca))
        {
            #-- row.sideColors
            if (!is.null(row.sideColors))
            {
                row.sideColors = as.matrix(row.sideColors)
                if (nrow(row.sideColors) != n)
                stop("'row.sideColors' must be a colors character vector
                (matrix) of length (nrow) ", n, ".",
                call. = FALSE)
            }
            
            #-- col.sideColors
            if (!is.null(col.sideColors))
            {
                col.sideColors = as.matrix(col.sideColors)
                if (nrow(col.sideColors) != p)
                stop("'col.sideColors' must be a colors character vector
                (matrix) of length (nrow) ", p, ".",
                call. = FALSE)
            }
            sample.sideColors = row.sideColors
            
            ## ----- DA object
            if(any(class.object %in%  c("mixo_splsda","mixo_plsda",
            'mixo_mlsplsda')))
            {
                #-- row.names
                if (is.logical(row.names))
                {
                    if (isTRUE(row.names))
                    row.names = mat$names$sample
                }
                #-- col.names
                if (is.logical(col.names))
                {
                    if (isTRUE(col.names))
                    col.names = mat$names$colnames$X
                }
                if(any(class.object %in%  c("mixo_splsda",'mixo_mlsplsda')))
                keep.X = apply(abs(mat$loadings$X[,comp, drop = FALSE]), 1,
                    sum) > 0
                else
                keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
                cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp],
                    use = "pairwise")
                X.mat = as.matrix(mat$variates$X[, comp])
            }
            else{
                #-- row.names
                if (is.logical(row.names))
                {
                    if (isTRUE(row.names))
                    row.names = mat$names$sample
                }
                #-- col.names
                if (is.logical(col.names))
                {
                    if (isTRUE(col.names))
                    col.names = mat$names$X
                }
                if(any(class.object %in%  c("spca","sipca")))
                keep.X = apply(abs(mat$rotation[,comp]), 1, sum) > 0
                else
                keep.X = apply(abs(mat$rotation), 1, sum) > 0
                cord.X = cor(mat$X[, keep.X], mat$x[, comp], use = "pairwise")
                X.mat = as.matrix(mat$x[, comp])
            }
            
            
            #-- cheking center and scale
            if (!is.logical(center)) {
                if (!is.numeric(center) || (length(center) != ncol(mat$X)))
                stop("'center' should be either a logical value or a numeric
                vector of length equal to the number of columns of 'X'.",
                call. = FALSE)
            }
            if (!is.logical(scale)) {
                if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
                stop("'scale' should be either a logical value or a numeric
                vector of length equal to the number of columns of 'X'.",
                call. = FALSE)
            }
            ## ---- clustering ----------------------------------------------
            object = scale(mat$X[, keep.X], center = center, scale = scale)
            col.names = col.names[keep.X]
            
            if (!is.null(col.sideColors))
            col.sideColors = as.matrix(col.sideColors[keep.X, ])
            
            if ((cluster == "both") || (cluster == "row")) {
                Rowv = rowMeans(X.mat)
                
                if (dist.method[1] == "correlation")
                dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)),
                method = "pearson"))
                else
                dist.mat = dist(X.mat, method = dist.method[1])
                
                hcr = hclust(dist.mat, method = clust.method[1])
                ddr = as.dendrogram(hcr)
                ddr = reorder(ddr, Rowv)
                rowInd = order.dendrogram(ddr)
                object = object[rowInd, ]
                row.names = row.names[rowInd]
                
                if (!is.null(row.sideColors))
                row.sideColors = as.matrix(row.sideColors[rowInd, ])
            }
            
            if ((cluster == "both") || (cluster == "column")) {
                Colv = rowMeans(cord.X)
                
                if (dist.method[2] == "correlation")
                dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                else
                dist.mat = dist(cord.X, method = dist.method[2])
                
                hcc = hclust(dist.mat, method = clust.method[2])
                ddc = as.dendrogram(hcc)
                ddc = reorder(ddc, Colv)
                colInd = order.dendrogram(ddc)
                object = object[, colInd]
                col.names = col.names[colInd]
                
                if (!is.null(col.sideColors))
                col.sideColors = as.matrix(col.sideColors[colInd, ])
            }
            
            #-- calling the image.map function --------------------------------#
            #------------------------------------------------------------------#
            
            
            
            #-- output --------------------------------------------------------#
            #------------------------------------------------------------------#
            res = list(mat = object, row.names = row.names,
            col.names = col.names)
            
            if ((cluster == "both") || (cluster == "row")) {
                res$rowInd = rowInd
                res$ddr = ddr
            }
            
            if ((cluster == "both") || (cluster == "column")) {
                res$colInd = colInd
                res$ddc = ddc
            }
            
            class(res) = paste("cim",class.object[1],sep="_")
            
        }
        else if(any(class.object %in%  object.rcc))
        {
            
            bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
            cord.X = cor(mat$X, bisect, use = "pairwise")
            cord.Y = cor(mat$Y, bisect, use = "pairwise")
            XY.mat = as.matrix(cord.X %*% t(cord.Y))
            
            #-- if mapping = "XY"
            if (mapping == "XY") {
                object = XY.mat
                
                cut=list()
                if (threshold != 0) {
                    cut[[1]] = unlist(lapply(1:nrow(object),
                    function(x){any(abs(object[x,]) > threshold)}))
                    object = object[cut[[1]],]
                    if (dist.method[1] != "correlation")
                    cord.X = cord.X[cut[[1]],]
                    
                    
                    if (is.null(nrow(object)) || nrow(object) == 0)
                    stop("threshold value very high. No variable was selected.",
                    call. = FALSE)
                    
                    
                    cut[[2]] = unlist(lapply(1:ncol(object),
                    function(x){any(abs(object[,x]) > threshold)}))
                    object = object[,cut[[2]]]
                    if (dist.method[2] != "correlation")
                    cord.Y = cord.Y[cut[[2]],]
                    
                    
                    if (is.null(ncol(object)) || ncol(object) == 0)
                    stop("threshold value very high. No variable was selected.",
                    call. = FALSE)
                    
                }
                
                if ((cluster == "both") || (cluster == "row")) {
                    #Rowv = rowMeans(XY.mat)
                    Rowv = rowMeans(cord.X)
                    
                    if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(object)),
                    method = "pearson"))
                    
                    else
                    dist.mat = dist(cord.X, method = dist.method[1])
                    
                    if (threshold > 0 ) {
                        row.names = row.names[cut[[1]]]
                        if (!is.null(row.sideColors))
                        row.sideColors = as.matrix(row.sideColors[cut[[1]], ])
                    }
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd, ]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd, ])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.Y)
                    
                    if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(as.matrix(object),
                    method = "pearson"))
                    else
                    dist.mat = dist(cord.Y, method = dist.method[2])
                    
                    if (threshold > 0 ) {
                        col.names = col.names[cut[[2]]]
                        if (!is.null(col.sideColors))
                        col.sideColors = as.matrix(col.sideColors[cut[[2]], ])
                    }
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd, ])
                }
            }
            
            #-- if mapping = "X"
            if (mapping == "X") {
                
                #-- cheking center and scale
                if (!is.logical(center)) {
                    if (!is.numeric(center) || (length(center) != ncol(mat$X)))
                    stop("'center' should be either a logical value or a
                    numeric vector of length equal to the number of columns of
                    'X'.",
                    call. = FALSE)
                }
                if (!is.logical(scale)) {
                    if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
                    stop("'scale' should be either a logical value or a
                    numeric vector of length equal to the number of columns of
                    'X'.",
                    call. = FALSE)
                }
                
                object = scale(mat$X, center = center, scale = scale)
                X.mat = as.matrix(mat$variates$X[, comp])
                
                if ((cluster == "both") || (cluster == "row")) {
                    Rowv = rowMeans(X.mat)
                    
                    if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)),
                    method = "pearson"))
                    else
                    dist.mat = dist(X.mat, method = dist.method[1])
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd, ]
                    row.names = row.names[rowInd]
                    
                    
                    if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd, ])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.X)
                    
                    if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                    else
                    dist.mat = dist(cord.X, method = dist.method[2])
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd, ])
                    
                }
                
            }
            
            #-- if mapping = "Y"
            if (mapping == "Y") {
                
                #-- cheking center and scale
                if (!is.logical(center)) {
                    if (!is.numeric(center) || (length(center) != ncol(mat$Y)))
                    stop("'center' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                    call. = FALSE)
                }
                if (!is.logical(scale)) {
                    if (!is.numeric(scale) || (length(scale) != ncol(mat$Y)))
                    stop("'scale' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                    call. = FALSE)
                }
                
                object = scale(mat$Y, center = center, scale = scale)
                Y.mat = as.matrix(mat$variates$Y[, comp])
                
                
                if ((cluster == "both") || (cluster == "row")) {
                    Rowv = rowMeans(Y.mat)
                    
                    if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(Y.mat)),
                    method = "pearson"))
                    else
                    dist.mat = dist(Y.mat, method = dist.method[1])
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd, ]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd, ])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.Y)
                    
                    if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(t(cord.Y), method = "pearson"))
                    else
                    dist.mat = dist(cord.Y, method = dist.method[2])
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd, ])
                    
                }
            }
            
            #-- calling the image.map function --------------------------------#
            #------------------------------------------------------------------#
            
            
            
            #-- output --------------------------------------------------------#
            #------------------------------------------------------------------#
            res = list(mat = object, row.names = row.names,
            col.names = col.names)
            
            if ((cluster == "both") || (cluster == "row")) {
                res$rowInd = rowInd
                res$ddr = ddr
            }
            
            if ((cluster == "both") || (cluster == "column")) {
                res$colInd = colInd
                res$ddc = ddc
            }
            
            class(res) = "cim_rcc"
            
        }
        else if(any(class.object %in%  object.pls))
        {
            if(any(class.object %in% c("mixo_spls","mixo_mlspls")))
            {
                keep.X = apply(abs(mat$loadings$X[,comp, drop = FALSE]), 1,
                sum) > 0
                keep.Y = apply(abs(mat$loadings$Y[,comp, drop = FALSE]), 1,
                sum) > 0}
            else
            {
                keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
                keep.Y = apply(abs(mat$loadings$Y), 1, sum) > 0}
            
            
            
            if (mat$mode == "canonical") {
                bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
                cord.X = cor(mat$X[, keep.X, drop = FALSE],
                             bisect, use = "pairwise")
                cord.Y = cor(mat$Y[, keep.Y, drop = FALSE],
                             bisect, use = "pairwise")
            }
            else {
                cord.X = cor(mat$X[, keep.X, drop = FALSE],
                mat$variates$X[, comp], use = "pairwise")
                cord.Y = cor(mat$Y[, keep.Y, drop = FALSE],
                mat$variates$X[, comp], use = "pairwise")
            }
            
            XY.mat = as.matrix(cord.X %*% t(cord.Y))
            sample.sideColors = row.sideColors
            
            #-- if mapping = "XY"
            if (mapping == "XY") {
                object = XY.mat
                row.names = row.names[keep.X]
                col.names = col.names[keep.Y]
                
                cut=list()
                if (threshold != 0) {
                    cut[[1]] = unlist(lapply(1:nrow(object),
                    function(x){any(abs(object[x,]) > threshold)}))
                    object = object[cut[[1]],]
                    if (dist.method[1] != "correlation")
                    cord.X = cord.X[cut[[1]],]
                    
                    
                    if (is.null(nrow(object)) || nrow(object) == 0)
                    stop("threshold value very high. No variable was selected.",
                    call. = FALSE)
                    
                    
                    cut[[2]] = unlist(lapply(1:ncol(object),
                    function(x){any(abs(object[,x]) > threshold)}))
                    object = object[,cut[[2]]]
                    if (dist.method[2] != "correlation")
                    cord.Y = cord.Y[cut[[2]],]
                    
                    
                    if (is.null(ncol(object)) || ncol(object) == 0)
                    stop("threshold value very high. No variable was selected.",
                    call. = FALSE)
                    
                }
                
                
                if (!is.null(row.sideColors))
                row.sideColors = as.matrix(row.sideColors[keep.X, ])
                if (!is.null(col.sideColors))
                col.sideColors = as.matrix(col.sideColors[keep.Y, ])
                
                if ((cluster == "both") || (cluster == "row")) {
                    
                    
                    #Rowv = rowMeans(XY.mat)
                    Rowv = rowMeans(cord.X)
                    
                    if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(object)),
                    method = "pearson"))
                    else
                    #dist.mat = dist(mat, method = dist.method[1])
                    dist.mat = dist(cord.X, method = dist.method[1])
                    
                    if (threshold > 0 ) {
                        row.names = row.names[cut[[1]]]
                        if (!is.null(row.sideColors))
                        row.sideColors = as.matrix(row.sideColors[cut[[1]], ])
                    }
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd, ]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd, ])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    #Colv = colMeans(mat)
                    Colv = rowMeans(cord.Y)
                    
                    if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(as.matrix(object),
                    method = "pearson"))
                    else
                    #dist.mat = dist(t(mat), method = dist.method[2])
                    dist.mat = dist(cord.Y, method = dist.method[2])
                    
                    if (threshold > 0 ) {
                        col.names = col.names[cut[[2]]]
                        if (!is.null(col.sideColors))
                        col.sideColors = as.matrix(col.sideColors[cut[[2]], ])
                    }
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd, ])
                }
            }
            
            #-- if mapping = "X"
            if (mapping == "X") {
                
                #-- cheking center and scale
                if (!is.logical(center)) {
                    if (!is.numeric(center) || (length(center) != ncol(mat$X)))
                    stop("'center' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'X'.",
                    call. = FALSE)
                }
                if (!is.logical(scale)) {
                    if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
                    stop("'scale' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'X'.",
                    call. = FALSE)
                }
                
                object = scale(mat$X[, keep.X], center = center, scale = scale)
                X.mat = as.matrix(mat$variates$X[, comp])
                col.names = col.names[keep.X]
                
                
                if (!is.null(col.sideColors))
                col.sideColors = as.matrix(col.sideColors[keep.X, ])
                
                
                
                if ((cluster == "both") || (cluster == "row")) {
                    Rowv = rowMeans(X.mat)
                    
                    if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)),
                    method = "pearson"))
                    else
                    dist.mat = dist(X.mat, method = dist.method[1])
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd, ]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd, ])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.X)
                    
                    if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                    else
                    dist.mat = dist(cord.X, method = dist.method[2])
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd, ])
                    
                }
            }
            
            #-- if mapping = "Y"
            if (mapping == "Y") {
                
                #-- cheking center and scale
                if (!is.logical(center)) {
                    if (!is.numeric(center) || (length(center) != ncol(mat$Y)))
                    stop("'center' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                    call. = FALSE)
                }
                if (!is.logical(scale)) {
                    if (!is.numeric(scale) || (length(scale) != ncol(mat$Y)))
                    stop("'scale' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                    call. = FALSE)
                }
                
                object = scale(mat$Y[, keep.Y], center = center, scale = scale)
                Y.mat = as.matrix(mat$variates$Y[, comp])
                col.names = col.names[keep.Y]
                
                
                if (!is.null(col.sideColors))
                col.sideColors = as.matrix(col.sideColors[keep.Y, ])
                
                
                
                if ((cluster == "both") || (cluster == "row")) {
                    Rowv = rowMeans(Y.mat)
                    
                    if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(Y.mat)),
                    method = "pearson"))
                    else
                    dist.mat = dist(Y.mat, method = dist.method[1])
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd, ]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd, ])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.Y)
                    
                    if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(t(cord.Y), method = "pearson"))
                    else
                    dist.mat = dist(cord.Y, method = dist.method[2])
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd, ])
                }
            }
            
            #-- calling the image.map function --------------------------------#
            #------------------------------------------------------------------#
            
            
            #-- output --------------------------------------------------------#
            #------------------------------------------------------------------#
            res = list(mat = object, row.names = row.names,
            col.names = col.names)
            
            if (!is.null(sample.sideColors)) {
                res$sample.sideColors = sample.sideColors
            }
            
            if ((cluster == "both") || (cluster == "row")) {
                res$rowInd = rowInd
                res$ddr = ddr
            }
            
            if ((cluster == "both") || (cluster == "column")) {
                res$colInd = colInd
                res$ddc = ddc
            }
            class(res) = paste("cim",class.object[1],sep="_")
        }} else {
    #-- if matrix class  -------------------------------------
            
        #-- mat
        isMat = tryCatch(is.matrix(mat), error = function(e) e)
        
        if ("simpleError" %in% class(isMat))
        stop(isMat[[1]], ".", call. = FALSE)
        
        if (!is.matrix(mat) || !is.numeric(mat))
        stop("'mat' must be a numeric matrix.", call. = FALSE)
        
        p = nrow(mat)
        q = ncol(mat)
        #-- row.names
        if (is.logical(row.names)) {
            if(isTRUE(row.names)) row.names = rownames(mat) else
            row.names = rep("", p)
        }
        else {
            row.names = as.vector(row.names)
            if (length(row.names) != p)
            stop("'row.names' must be a character vector of length ", p, ".",
            call. = FALSE)
        }
        
        #-- col.names
        if (is.logical(col.names)) {
            if(isTRUE(col.names)) col.names = colnames(mat) else
            col.names = rep("", q)
        }
        else {
            col.names = as.vector(col.names)
            if (length(col.names) != q)
            stop("'col.names' must be a character vector of length ", q, ".",
            call. = FALSE)
        }
        
        #-- row.sideColors
        if (!is.null(row.sideColors)) {
            row.sideColors = as.matrix(row.sideColors)
            if (nrow(row.sideColors) != p)
            stop("'row.sideColors' must be a colors character vector (matrix)
            of length (nrow) ", p, ".",
            call. = FALSE)
        }
        
        #-- col.sideColors
        if (!is.null(col.sideColors)) {
            col.sideColors = as.matrix(col.sideColors)
            if (nrow(col.sideColors) != q)
            stop("'col.sideColors' must be a colors character vector (matrix)
            of length (nrow) ", q, ".",
            call. = FALSE)
        }
        
        #-- clustering -------------------------------------------------------#
        #---------------------------------------------------------------------#
        object=mat
        if ((cluster == "both") || (cluster == "row")) {
            Rowv = rowMeans(mat)
            
            if (dist.method[1] == "correlation")
            dist.mat = as.dist(1 - cor(t(as.matrix(mat)), method = "pearson"))
            else
            dist.mat = dist(mat, method = dist.method[1])
            
            hcr = hclust(dist.mat, method = clust.method[1])
            ddr = as.dendrogram(hcr)
            ddr = reorder(ddr, Rowv)
            rowInd = order.dendrogram(ddr)
            object = mat[rowInd, ]
            row.names = row.names[rowInd]
            
            if (!is.null(row.sideColors))
            row.sideColors = as.matrix(row.sideColors[rowInd, ])
        }
        
        if ((cluster == "both") || (cluster == "column")) {
            Colv = colMeans(mat)
            
            if (dist.method[2] == "correlation")
            dist.mat = as.dist(1 - cor(as.matrix(mat), method = "pearson"))
            else
            dist.mat = dist(t(mat), method = dist.method[2])
            hcc = hclust(dist.mat, method = clust.method[2])
            ddc = as.dendrogram(hcc)
            ddc = reorder(ddc, Colv)
            colInd = order.dendrogram(ddc)
            object = object[, colInd]
            col.names = col.names[colInd]
            
            if (!is.null(col.sideColors))
            col.sideColors = as.matrix(col.sideColors[colInd, ])
        }
        #-- calling the image.map function ------------------------------------#
        
        
        
        #-- output ------------------------------------------------------------#
        #----------------------------------------------------------------------#
        res = list(mat = object, row.names = row.names, col.names = col.names,
        row.sideColors = row.sideColors, col.sideColors = col.sideColors)
        
        if ((cluster == "both") || (cluster == "row")) {
            res$rowInd = rowInd
            res$ddr = ddr
        }
        
        if ((cluster == "both") || (cluster == "column")) {
            res$colInd = colInd
            res$ddc = ddc
        }
        
        class(res) = "cim_default"
        
    }
    #-- call imageMap  -------------------------------------
    opar = par(no.readonly = TRUE)
    
    try_plot <- tryCatch({
    imageMap(object,
    color = color,
    row.names = row.names,
    col.names = col.names,
    row.sideColors = row.sideColors,
    col.sideColors = col.sideColors,
    row.cex = row.cex,
    col.cex = col.cex,
    cluster = cluster,
    ddr = ddr,
    ddc = ddc,
    cut.tree = cut.tree,
    transpose = transpose,
    symkey = symkey,
    keysize = keysize,
    keysize.label = keysize.label,
    zoom = zoom,
    title = title,
    xlab = xlab,
    ylab = ylab,
    margins = margins,
    lhei = lhei,
    lwid = lwid)}, error = function(e) e)

    if (is(try_plot, "error")) {
        message(sprintf("Error in cim: %s. Try expanding the plot pane in RStudio and deleteing existing plots. If this error persists put save=TRUE, name.save = 'a_name_for_plot.pdf' to save the plot.", try_plot$message))
    } else {
        #-- add to plot  -------------------------------------
        if (!is.null(legend))
        {if(is.null(legend$x)) legend$x = "topright"
        if(is.null(legend$bty)) legend$bty = "n"
        if (is.null(legend$cex)) legend$cex = 0.8
        if(any(class.object %in%  c("mixo_splsda","mixo_plsda")))
        {
            if (is.null(legend$legend)) legend$legend = mat$names$colnames$Y
            
            #-- col
            if (is.null(legend$col)) {
                if (!is.null(sample.sideColors))
                    legend$col = unique(as.matrix(sample.sideColors[order(
                        map(mat$ind.mat)), 1]))
            }
            
        }
        else if(any(class.object %in%  c("mixo_mlsplsda")))
        {
            if (is.null(legend$legend) && is.null(legend$col)) {
                if (ncol(mat$multilevel) >= 2) {
                    df = data.frame(mat$multilevel[, 2], sample.sideColors[, 1])
                    df = unique(df)
                    legend$legend = as.character(df[, 1])
                    legend$col = as.character(df[, 2])
                }
                if (ncol(mat$multilevel) == 3) {
                    df = data.frame(mat$multilevel[, 3], sample.sideColors[, 2])
                    df = unique(df)
                    legend$legend = c(legend$legend, as.character(df[, 1]))
                    legend$col = c(legend$col, as.character(df[, 2]))
                }
            }
        }
        else if(any(class.object %in%  c("mixo_mlspls")))
        {
            if (mapping != "XY") {
                if (is.null(legend$legend) && is.null(legend$col)) {
                    if (ncol(mat$multilevel) >= 2) {
                        df = data.frame(mat$multilevel[, 2],
                                        sample.sideColors[, 1])
                        df = unique(df)
                        legend$legend = as.character(df[, 1])
                        legend$col = as.character(df[, 2])
                    }
                    if (ncol(mat$multilevel) == 3) {
                        df = data.frame(mat$multilevel[, 3],
                                        sample.sideColors[, 2])
                        df = unique(df)
                        legend$legend = c(legend$legend, as.character(df[, 1]))
                        legend$col = c(legend$col, as.character(df[, 2]))
                    }
                }
            }
        }
        if (is.null(legend$legend))
            stop("argument \"legend$legend\" is missing, with no default")
        
        #-- fill
        if (is.null(legend$fill)) legend$fill = legend$col
        
        
        par(mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
        
        if (!is.null(legend$title))
        {
            legend(x = legend$x, y = legend$y, legend = legend$legend,
                   col = legend$col, fill = legend$fill, bty = legend$bty,
                   title=legend$title,cex=legend$cex)
        }else{
            legend(x = legend$x, y = legend$y, legend = legend$legend,
                   col = legend$col, fill = legend$fill, bty = legend$bty,
                   cex=legend$cex)
        }
        
        
        }
        if (any(class.object %in% object.list) & !any(class.object %in%
                                                      object.pca) & mapping =="XY")
            res$mat.cor=object
    }
    
    par(opar)
    
    if (!is.null(save))
    dev.off()
    
    return(invisible(res))
}

