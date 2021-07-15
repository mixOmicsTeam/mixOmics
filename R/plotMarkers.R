#' Plot the values for multivariate markers in block analyses
#'
#' Plots the standardised values (after centring and/or scaling) for the
#' selected variables for a given block on a given component. Only applies to
#' \code{block.splsda} or \code{block.spls}.
#' 
#' @param object An object of class \code{block.splsda} or \code{block.spls}
#' @param block Name or index of the block to use
#' @param comp Integer, the component to use
#' @param markers Character or integer, only include these markers. If integer,
#'   the top 'markers' features are shown
#' @param group Factor, the grouping variable (only required for
#'   \code{block.spls} objects)
#' @template arg/col.per.group
#' @param global Logical indicating whether to show the global plots (TRUE) or
#'   segregate by feature (FALSE)
#' @param title The plot title
#'
#' @return A ggplot object
#' @seealso \code{\link{plotLoadings}}, \code{\link{block.splsda}}, \code{\link{block.spls}}
#' @export
#'
#' @examples
#' # see ?block.splsda and ?block.spls
plotMarkers <-
    function(object,
             block,
             markers = NULL,
             comp = 1,
             group = NULL,
             col.per.group = NULL,
             global = FALSE,
             title = NULL)
    {
        
    blocks <- names(object$X)
    if (is.numeric(block))
    {
        blocks <- blocks[block]
    }
    if (!block %in% blocks)
        stop(message = sprintf("block must be one of: %s", paste0(blocks, collapse = ', ')))
    df <- data.frame(object$X[[block]], check.names = FALSE)
    
    ## group factor
    group <- .get.group(group, object, n_ind = nrow(df))
    
    col.group <- .get.cols.and.group(col.per.group = col.per.group, 
                                     group = group)
    group <- col.group$group
    col.per.group <- col.group$col.per.group
    vars <- selectVar(object, block=block, comp=comp)[[1]]$value
   
    var.names <- rownames(vars)
    
    df <- df[,var.names]
    df$group <- group
    df <-
        melt(df,
             id.vars = 'group',
             variable.name = 'feature',
             value.name = 'value')
    df$feature <- factor(df$feature, levels = var.names, ordered = TRUE)
    df$sign <- ifelse(df$value > 0, 'Positive Loading', 'Negative Loading')
    # to show +ves on top
    df$sign <- factor(df$sign, levels = c('Positive Loading', 'Negative Loading'), ordered = TRUE)
    
    if (!is.null(markers))
    {
        if (is.numeric(markers))
        {
            markers <- selectVar(object = object, comp = comp)[[block]]$name
        }
        if ( is.character(markers))
        {
            invalid.markers <- setdiff(markers, df$feature)
            if (length(invalid.markers) > 0)
                stop("invalid feature names: ", paste0(invalid.markers, collapse = ", "), call. = FALSE)
            #' @importFrom dplyr filter
            feature <- NULL
            df <- filter(df, feature %in% markers)
        }
    }
    if (global)
    {
        df$feature <- 'plotMarkers'
    }
    
    if (is.null(title))
        title <- sprintf("Block: %s | Component: %s", block, comp)
    p <- ggplot(df, aes_string('group', 'value', fill='group')) +  
        geom_violin(adjust=0.9) +
        geom_boxplot(width=0.1) + 
        scale_fill_manual(values = col.per.group) + 
        theme_classic() + 
        labs(x='', 
             y='value (standardised)', 
             title = title) + 
        theme(legend.position = 'none', 
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
              plot.title = element_text(hjust = 0.5, colour = "grey40"), 
              strip.background = element_rect(colour="black", fill="grey80"))
    if (global)
    {
        p <- p + facet_grid(sign~feature, scales = 'free') 
    } else {
        p <- p + facet_grid(.~feature, scales = 'free') 
    }
    p
}
