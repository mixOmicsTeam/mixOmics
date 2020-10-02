
## ------------------------------ plot.spca ------------------------------- ##
#' Show (s)pca explained variance plots
#' 
#' @param x A \code{(s)pca} object
#' @param ncomp Integer, the number of components
#' @param type Character, default "barplot" or any other type available in plot, as "l","b","p",..
#' @param ... Not used
#' @author Kim-Anh Lê Cao, Florian Rohart, Leigh Coonan, Al J Abadi
#' @method plot pca
#' @export
plot.pca  <- function(x,
                      ncomp = length(x$explained_variance),
                      type = "barplot",
                      # either barplot or any other type available in plot, as "l","b","p",..
                      ...)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- ncomp
    ncomp_model <- length(x$explained_variance)
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
        stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > ncomp_model)
        stop("'ncomp' must be lower or equal to ", length(ncomp_model), ".",
             call. = FALSE)
    ## end check - begin screeplot
    expl_vars = (x$explained_variance)[seq_len(ncomp)] # relative variance
    ylab = "Explained Variance"
    if (type == "barplot")
    {
        barplot(expl_vars, names.arg = seq_len(ncomp), xlab = "Principal Components", ylab = ylab,...)
    } else {
        plot(expl_vars, type = type, axes = FALSE,
             xlab = "Principal Components",
             ylab = ylab,... )
        axis(1, at = seq_len(ncomp))
        axis(2)
    }
    
}

## ------------------------------- plot.rcc ------------------------------- ##
#' Canonical Correlations Plot
#' 
#' This function provides scree plot of the canonical correlations.
#' 
#' @inheritParams plot.pca
#' @param x object of class inheriting from \code{"rcc"}.
#' @return none
#' @author Sébastien Déjean, Ignacio González, Al J Abadi
#' @seealso \code{\link{points}}, \code{\link{barplot}}, \code{\link{par}}.
#' @keywords multivariate hplot
#' @method plot rcc
#' @examples
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)
#' 
#' ## 'pointplot' type scree
#' plot(nutri.res) #(default)
#' 
#' \dontrun{
#' plot(nutri.res, pch = 19, cex = 1.2,
#' col = c(rep("red", 3), rep("darkblue", 18)))
#' 
#' ## 'barplot' type scree
#' plot(nutri.res, type = "barplot")
#' 
#' plot(nutri.res, type = "barplot", density = 20, col = "black")
#' }
#' 
#' @export
plot.rcc <-
    function(x, type = "barplot", ...) 
    {
        if (hasArg(scree.type)) {
           stop("'scree.type' has been replaced by 'type'. See ?plot.rcc\n")
        }
        if (type == "barplot") {
            barplot(x$cor, xlab = "Dimension", ylim = c(0, 1),
                    ylab = "Canonical correlation", ...)
        }
        else {
            plot(x$cor, xlab = "Dimension", ylim = c(0, 1),
                 ylab = "Canonical correlation", type = type, ...)
        }
    }
