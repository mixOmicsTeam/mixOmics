#' biplot methods for \code{pca} family
#'
#' @inheritParams plotIndiv
#' @inheritParams plotVar
#' @param x A 'pca' object.
# @param hide One of \code{c('none', 'ind', 'var')} indicating whether samples
#   or variables should not be shown
#' @param ind.names.size Numeric, sample name size.
#' @param ind.names.col Character, sample name colour.
#' @param ind.names.repel Logical, whether to repel away label names.
#' @param group Factor indicating the group membership for each sample.
#' @param pch.size Numeric, sample point character size.
#' @param pch.levels If \code{pch} is a factor, a named vector providing the
#'   point characters to use. See examples.
#' @param var.names Logical indicating whether to show variable names. 
#' Alternatively, a character.
#' @param var.names.col Character, variable name colour.
#' @param var.names.size Numeric, variable name size.
#' @param var.names.angle Logical, whether to align variable names to arrow
#'   directions.
#' @param var.arrow.col Character, variable arrow colour. If 'NULL', no arrows
#'   are shown.
#' @param var.arrow.size Numeric, variable arrow head size.
#' @param var.arrow.length Numeric, length of the arrow head in 'cm'.
#' @param ind.legend.title Character, title of the legend.
#' @param vline Logical, whether to draw the vertical neutral line.
#' @param hline Logical, whether to draw the horizontal neutral line.
# @param block Character, currently the only valid value is 'X'.
#' @param legend Logical, whether to show the legend if \code{group != NULL}.
#' @param legend.title Character, the legend title if \code{group != NULL}.
#' @param pch.legend Character, the legend title if \code{pch} is a factor.
#' @param pch.legend.title Character, the legend title if \code{pch} is a factor.
#' @param ... Not currently used.
#' @details 
#' \code{biplot} unifies the reduced representation of both the
#' observations/samples and variables of a matrix of multivariate data on the
#' same plot. Essentially, in the reduced space the samples are shown as
#' points/names and the contributions of features to each dimension are shown as
#' directed arrows or vectors.
#' @return A ggplot object.
#' @author Al J Abadi
#' @importFrom ggrepel geom_text_repel
#' @example ./examples/biplot-examples.R
#' @name biplot
#' @method biplot pca
#' @export
biplot.pca <- function(x,
                       comp = c(1,2),
                       # hide = c('none', 'ind', 'var'), ## should not be necessary 
                       ind.names = TRUE,
                       group = NULL,
                       cutoff = 0,
                       col.per.group=NULL,
                       col = NULL,
                       ind.names.size = 3,
                       ind.names.col = color.mixo(4),
                       ind.names.repel = TRUE,
                       pch = 19,
                       pch.levels=NULL,
                       pch.size = 2,
                       var.names = TRUE,
                       var.names.col = 'grey40',
                       var.names.size = 4,
                       var.names.angle = FALSE,
                       var.arrow.col = 'grey40',
                       var.arrow.size = 0.5,
                       var.arrow.length = 0.2,
                       ind.legend.title = NULL,
                       vline = FALSE,
                       hline = FALSE,
                       # block = 'X', # also for PLS methods?
                       legend = if (is.null(group)) FALSE else TRUE,
                       legend.title = NULL,
                       pch.legend.title = NULL,
                       ...
)
{
    object <- x
    rm(x)
    block <- 'X'
    # hide <- match.arg(hide)
    hide <- 'none'
    selection <- rowSums(object$loadings$X[, comp]) != 0 
    loadings <- object$loadings[[block]][selection, comp]
    loadings <- data.frame(loadings)
    
    ## scale check
    if (isFALSE(object$call$scale))
        warning("The 'tune.spca' algorithm has used scale=FALSE. We recommend scaling the data",
                " to improve orthogonality in the sparse components.")
    ## cutoff correlation
    cutoff <- .check.cutoff(cutoff)
    cors <- cor(object$X[, selection], object$variates$X[, comp], use = 'pairwise' )
    above.cutoff <- apply(cors, 1, function(x) any(abs(x) >= cutoff))
    loadings <- loadings[above.cutoff,]
    
    variates <- object$variates[[block]][, comp]
    variates <- data.frame(variates)
    ## scaler of var vs sample coordinates
    scaler <- max(variates, na.rm = TRUE)/max(loadings, na.rm = TRUE)
    
    PCs <- colnames(variates)
    expl_vars <- round(object$explained_variance*100)[comp]
    axes.titles <- sprintf("%s (%s%%)", PCs, expl_vars)
    ind.names <- .get.character.vector(ind.names, vec = rownames(variates))
    
    variates$ind.names <- ind.names
    col.group <-
        .get.cols.and.group(
            col.per.group = col.per.group,
            group = group,
            col = col,
            object = object,
            n_ind = nrow(variates)
        )
    group <- col.group$group
    col.per.group <- col.group$col.per.group
    if (length(col.per.group) == 1)
    {
        legend <- FALSE
    }
    
    ## ------------- outline
    gg_biplot <- 
        ggplot() + 
        theme_classic() +  
        labs(x = axes.titles[1], 
             y = axes.titles[2])
    ## vline and hline
    if (vline)
    {
        gg_biplot <- gg_biplot + geom_vline(xintercept = 0, size = 0.3, col = 'grey75')
    }
    if (hline)
    {
        gg_biplot <- gg_biplot +  geom_hline(yintercept = 0, size = 0.3, col = 'grey75')
    }
   
        
    ## ------------- inds
    if (! 'ind' %in% hide) 
    {
        if (!is.null(pch)) 
        {
            ## ------------- advanced user args
            fill <- ifelse(is.null(list(...)$fill), 'black', list(...)$fill)
            alpha <- ifelse(is.null(list(...)$alpha), 1, list(...)$alpha)
            
            pch.res <- .get.pch(pch, pch.levels, n_ind = nrow(variates))
            pch <- pch.res$pch
            pch.levels <- pch.res$pch.levels
            pch.legend <- pch.res$pch.legend
            
            ## get 'pch' and 'group' arg used for legends so we can handle
            ## legends whether needed or not in a unified way (see scale_*_manual)
            pch.legend.title <- .change_if_null(pch.legend.title, as.character(as.list(match.call())['pch']))
            legend.title <- .change_if_null(legend.title, as.character(as.list(match.call())['group']))
            
            # pch.col <- .get.ind.colors(group, col, col.per.group, n_ind = nrow(variates))
            gg_biplot <- gg_biplot + 
                geom_point(aes(x = variates[, comp[1]], 
                               y = variates[, comp[2]],
                               col = group,
                               shape = pch),
                           fill = fill,
                           alpha = alpha,
                           size = pch.size) + 
                scale_shape_manual(values = pch.levels,
                                   guide = if (isTRUE(pch.legend)) guide_legend(title = pch.legend.title) else NULL)
            
        }
        else
        {
            ind.names.repel <- FALSE
        }
        if (!is.null(ind.names))
        {
            if (isTRUE(ind.names.repel)) {
                gg_biplot <- gg_biplot + 
                    geom_text_repel(mapping = aes(x = variates[, comp[1]],
                                                  y = variates[, comp[2]],
                                                  label = ind.names,
                                                  col = group
                                    ), 
                                    size = ind.names.size,
                                    show.legend = FALSE)
            } else {
                gg_biplot <- gg_biplot + 
                    geom_text(mapping = aes(x = variates[, comp[1]],
                                                  y = variates[, comp[2]],
                                                  label = ind.names,
                                                  col = group
                    ), 
                    size = ind.names.size,
                    show.legend = FALSE)
            }
            
        }
        gg_biplot <- gg_biplot + 
            scale_color_manual(values = col.per.group,
                               guide = if (isTRUE(legend)) guide_legend(title = legend.title) else NULL)
    }
    
    ## ------------- vars
    
    if (! 'var' %in% hide) 
    {
        loadings <- loadings*scaler
        var.names.col <- .get.ind.colors(group = NULL, 
                                         col = var.names.col,
                                         col.per.group = NULL, 
                                         n_ind = nrow(loadings))
        if (!is.null(var.arrow.col))
        {
            var.arrow.col <- .get.ind.colors(group = NULL, 
                                             col = var.arrow.col,
                                             col.per.group = NULL, 
                                             n_ind = nrow(loadings))
            loadings$var.names.col <- var.names.col
            loadings$var.arrow.col <- var.arrow.col
            ## lines and arrows
            gg_biplot <-
                gg_biplot + geom_segment(
                    aes(
                        x = 0,
                        y = 0,
                        xend = loadings[,comp[1]],
                        yend = loadings[,comp[2]],
                    ),
                    col = var.arrow.col,
                    arrow = arrow(length = unit(var.arrow.length, "cm")),
                    size = var.arrow.size,
                    show.legend = FALSE
                )
        }
        
        ## labels
        var.labels <- .get.character.vector(arg = var.names, vec = rownames(loadings))
        ## label angles
        angle <- rep(0, nrow(loadings))
        
        if (!is.null(var.names)) 
        {
            angle <- rep(0, nrow(loadings))
            if (var.names.angle == TRUE)
            {
                angle <- atan(loadings[, comp[2]]/loadings[, comp[1]]) * 360/(2 * pi)
            }
            
            gg_biplot <-
                gg_biplot + geom_text_repel(
                    aes(
                        x = loadings[, comp[1]],
                        y = loadings[, comp[2]],
                        label = var.labels,
                        angle = angle,
                        hjust = ifelse(loadings[, comp[1]] > 0, 1, 0),
                        vjust = ifelse(loadings[, comp[2]] > 0, 1, 0)
                    ),
                    col = var.names.col,
                    size = var.names.size,
                    box.padding = 0.1,
                    point.padding = 0.1
                )
        } 
        
        ## second set of axes
        gg_biplot <- gg_biplot + scale_y_continuous(sec.axis = sec_axis(~.*1/scaler)) +
            scale_x_continuous(sec.axis = sec_axis(~.*1/scaler)) 
    }
    gg_biplot
}
