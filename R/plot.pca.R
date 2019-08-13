#' Plot the explained variances from a pca object
#'
#' Creates a scree plot of explained variance by the study PCs.
#'
#' @title plot methods for mixOmics

## ----------------------------------- Parameters
#' @param x a \code{pca} object obtained from \code{pca} function.
#' @param ncomp number of PCs to show.
#' @param type type of the plot, either "barplot" or argument passed to \code{type} in base \code{plot}.
#' @param explained.var logical. Whether to show proportion of variance explained (TRUE) or the total variance (FALSE).
#' @param ... other arguments passed to \code{plot}.

## ----------------------------------- Value
#' @return A scree plot of explained variance by the study PCs.
## ----------------------------------- Misc
#' @author Florian Rohart, Kim-Anh Lê Cao, Ignacio González, Al J Abadi
#' \code{\link{plotIndiv}}, \code{\link{pca}} and http://www.mixOmics.org
#' for more details.
#' @keywords plot

## ----------------------------------- Examples
#' @example examples/plot-example.R

## ----------------------------------- Method
#' @rdname plot
#' @export
plot.pca <-  function(x,
            ncomp = min(10, length(x$sdev)),
            type = "barplot", # either barplot or any other type available in plot, as "l","b","p",..
            explained.var=TRUE,
            ...)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#

    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)

    ncomp = round(ncomp)

    if (ncomp > length(x$sdev))
    stop("'ncomp' must be lower or equal than ", length(x$sdev), ".",
    call. = FALSE)

    #-- end checking --#
    #------------------#

    #-- scree plot -------------------------------------------------------------#
    #---------------------------------------------------------------------------#

    variances = (x$sdev^2)[1:ncomp] # relative variance
    ylab = "Variance"
    if(explained.var==TRUE)
    {
        variances=variances/x$var.tot #explained variances
        ylab = "Explained Variance"
    }
    if (type == "barplot")
    {
        barplot(variances, names.arg = seq(1, ncomp), xlab = "Principal Components", ylab = ylab,...)
    } else {
        plot(variances, type = type, axes = FALSE,
        xlab = "Principal Components",
        ylab = ylab,... )
        axis(1, at = 1:ncomp)
        axis(2)
    }

}


