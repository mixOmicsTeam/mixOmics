## ------------------------ .plotLoadings_barplot ------------------------- ##
#' plotLoadings helper
#'
#' @noRd
.plotLoadings_barplot <- function(height, col, names.arg, cex.name, border, xlim) {
    tryCatch({barplot(height, horiz = TRUE, las = 1, col = col, axisnames = TRUE, names.arg = names.arg, #names.arg = row.names(df),
                      cex.names = cex.name, cex.axis = 0.7, beside = TRUE, border = border, xlim = xlim)},
             error = function(e){
                 if ( grepl(pattern = "figure margins too large", e) ){
                     stop("\nplotLoadings encountered margin errors. Ensure feature names are not too long (see 'name.var' argument) and the 'Plots' pane is cleared and enlargened.\n", call. = FALSE)
                 } else {
                     stop(e$message, call. = FALSE)
                 }
             })
}

## ----------------------------- .check_cpus ------------------------------ ##
#' Check cpus argument
#'
#' @param cpus Integer, the number of cpus for parallel processing.
#'
#' @return The number of cpus available/applicable, or condition.
#' @author Al J Abadi
#'
#' @noRd
.check_cpus <- function(cpus) {
    if (!is.numeric(cpus) || cpus <= 0)
        stop("Number of CPUs should be a positive integer.", call. = FALSE)
    
    if (cpus < 0 || is.null(cpus))
        cpus <- 1
    
    if (cpus > parallel::detectCores())
        message(sprintf("\nOnly %s CPUs available for parallel processing.\n", parallel::detectCores()))
    return(cpus)
    
}

## ----------------------------- .check_alpha ----------------------------- ##
#' Check significance threshold sanity
#'
#' @param alpha numeric, significance threshold for t.test
#'
#' @return NULL if acceptable value, otherwise condition
#' @author Al J Abadi
#'
#' @noRd
.check_alpha <- function(alpha=NULL) {
    if (is.null(alpha))
        alpha <- 0.01
    
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1)
        stop("invalid 'signif.threshold'. Use 0.01 or 0.05.", call. = FALSE)
    alpha
}

## --------------------------- .unexpected_err ---------------------------- ##
#' Unexpected error handler for the package
#'
#' To be used in unexperimented situations where handlers fail
#' @param trying_to character, the context in which unexpected error occurred.
#'
#' @return Error message
#' @author Al J Abadi
#'
#' @noRd
.unexpected_err <- function(trying_to = NULL) {
    trying_to <- ifelse(is.null(trying_to), "", sprintf(" while trying to %s", trying_to))
    msg <- sprintf("Unexpected error%s. Please check the inputs and if problem persists submit an issue to https://github.com/mixOmicsTeam/mixOmics/issues", trying_to)
    stop(msg, call. = FALSE)
}

## ------------------------------- .on_unix ------------------------------- ##
#' Check OS type for parallel processing
#'
#' @return Logical, FALSE if windows OS, TRUE if unix OS.
#' @author Al J Abadi
#'
#' @noRd
.onUnix <- function() {
    return(ifelse(.Platform$OS.type == "unix", TRUE, FALSE))
}

## ----------------------- .unlist_repeat_cv_output ----------------------- ##
#' repeat_cv_perf.diablo helper to unlist internal outputs to
#' previous list format for downstream work
#'
#' @return list of outputs needed for downstream work
#' @author Al J Abadi
#'
#' @noRd
.unlist_repeat_cv_output <- function(list_rep=NULL) {
    
    ## so lapply will return a named list
    names_list <- as.list(names(list_rep[[1]]))
    names(names_list) <- names(list_rep[[1]])
    
    list_nrep <- lapply(names_list, function(x){
        lapply(list_rep, function(y) {y[[x]]})
    })
    
    
    return(list_nrep)
}

## ------------------------ stratified_subsampling ------------------------ ##
#' Perform stratified subsampling for cross-validation
#'
#' @param Y A factor or a class vector for the discrete outcome
#' @param folds  Integer, if validation = Mfold, how many folds?
#' Mainly to be used for reproducibility and unit tests.
#'
#' @return A list containing:
#' \item{SAMPLE} List of samples for each fold
#' \item{stop} A warning signal for when there is not enough samples in a 
#' class
#' @author Florian Rohart, Al J Abadi
#' @noRd
stratified.subsampling <- function(Y, folds = 10)
{
    
    stop <- 0
    for (i in seq_len(nlevels(Y)))
    {
        ai <- sample(which(Y == levels(Y)[i]), replace = FALSE) # random sampling of the samples from level i
        aai <- suppressWarnings(split(ai, factor(seq_len(
            min(folds, length(ai))
        ))))                       # split of the samples in k-folds
        if (length(ai) < folds)
            # if one level doesn't have at least k samples, the list is completed with "integer(0)"
        {
            for (j in (length(ai) + 1):folds)
                aai[[j]] <- integer(0)
            stop <- stop + 1
        }
        assign(paste("aa", i, sep = "_"), sample(aai, replace = FALSE))         # the `sample(aai)' is to avoid the first group to have a lot more data than the rest
    }
    
    # combination of the different split aa_i into SAMPLE
    SAMPLE <- list()
    for (j in seq_len(folds))
    {
        SAMPLE[[j]] <- integer(0)
        for (i in seq_len(nlevels(Y)))
        {
            SAMPLE[[j]] <- c(SAMPLE[[j]], get(paste("aa", i, sep = "_"))[[j]])
        }
    }# SAMPLE is a list of k splits
    
    ind0 <- sapply(SAMPLE, length)
    if (any(ind0 == 0))
    {
        SAMPLE <- SAMPLE [-which(ind0 == 0)]
        message(
            "Because of a too high number of 'folds' required, ",
            length(which(ind0 == 0)),
            " folds were randomly assigned no data: the number of 'folds' is reduced to ",
            length(SAMPLE)
        )
    }
    
    return(list(SAMPLE = SAMPLE, stop = stop))
}

## ------------------------------ .name_list ------------------------------ ##
#' Create a named list
#'
#' Creates a named list of a character vector where names and values
#' are the same, useful for lapply family.
#' @param char Character vector.
#'
#' @return A named list.
#' @author Al J Abadi
#'
#' @noRd
.name_list <- function(char, names=NULL) {
    out <- as.list(char)
    names <- if (is.null(names)) char else names
    names(out) <- names
    return(out)
}

## ------------------------------ .mixo_rng ------------------------------- ##
#' RNG used for package tests
#'
#' Creating a function so that it can easily be changed.
#' It's mainly used in form of the given example in unit tests.
#'
#' @return An R version
#' @author Al J Abadi
#' @example RNGversion(.mixo_rng())
#' @noRd
.mixo_rng <- function() {"3.6.0"}

## --------------------------- .change_if_null ---------------------------- ##
#' Set default value if NULL
#'
#' For arguments with default NULL value, sets the desired default value
#' 
#' @param arg The function argument with NULL/no default value
#' @param default The desired default value
#'
#' @return The desired default value
#' @author Al J Abadi
#'
#' @examples
#' .change_if_null(NULL, 10)
#' #> 10
#' .change_if_null(3, 10)
#' #> 3
#' @noRd
.change_if_null <- function(arg, default) {
    return(
        if ( missing(arg) || is.null(arg) )
            default
        else
            arg
    )
}

## ------------------------ .add_consensus_blocks ------------------------- ##
#' Add consensus blocks to DIABLO object
#'
#' For plotIndiv(..., blocks = "consensus")
#' 
#' @param block_object A diablo object
#' @param consensus_blocks One or both of c('consensus', 'weighted.consensus')
#'
#' @return A diablo object with consensus blocks added for smooth input into plotIndiv.sgccda
#' @author Al J Abadi
#' @examples 
#' \dontrun{
#' data(nutrimouse)
#' Y = nutrimouse$diet
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
#' design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)
#' 
#' nutrimouse.sgccda1 <- wrapper.sgccda(X = data,
#'                                      Y = Y,
#'                                      design = design1,
#'                                      ncomp = 2,
#'                                      keepX = list(gene = c(10,10), lipid = c(15,15)),
#'                                      scheme = "centroid")
#' nutrimouse.sgccda.consensus <- mixOmics:::.add_consensus_blocks(nutrimouse.sgccda1, consensus_blocks = "consensus")
#' names(nutrimouse.sgccda.consensus$X)
#> "gene"      "lipid"     "consensus"
#' }
#' @noRd
.add_consensus_blocks <- function(block_object, consensus_blocks = c('consensus', 'weighted.consensus')) {
    X_blocks <- with(block_object, names$blocks[-which(names$block == 'Y')])
    consensus_blocks <- match.arg(consensus_blocks, several.ok = TRUE)
    
    .get_consensus_variates <- function(object, X_blocks, weighted = FALSE) {
        
        arrays <- object$variates
        arrays <- arrays[X_blocks]
        if (isTRUE(weighted)) {
            if (!is(block_object, "sgccda")) {
                if ('weighted.consensus' %in% consensus_blocks ) {
                    stop("'weighted.consensus' plots are only available for block.splsda objects ")
                }
            }
            weights <- object$weights
        } else {
            ## if unweighted, create a weight data.frame of 1 to avoid code complications
            weights <- matrix(rep(1, length(X_blocks)*object$ncomp[1]), nrow = length(X_blocks))
            dimnames(weights) <- list(X_blocks, paste0("comp", seq_len(object$ncomp[1])))
            weights <- as.data.frame(weights)
        }
        block_names <- .name_list(names(arrays))
        
        weighted_arrays <- lapply(block_names, function(x){
            variates <- arrays[[x]]
            weights <- diag(weights[x, ])
            weighted_variates <- variates %*% weights
            dimnames(weighted_variates) <- dimnames(variates)
            weighted_variates
        })
        wtd_sum <- Reduce(f = '+', weighted_arrays)
        ## weighted mean = weighted sum / sum(weight)
        sweep(wtd_sum, MARGIN = 2, colSums(weights), FUN = "/")
    }
    
    for (consensus_block in consensus_blocks) {
        block_object$X[[consensus_block]] <-  0
        if (consensus_block == "weighted.consensus") {
            block_object$variates[[consensus_block]] <-  .get_consensus_variates(object = block_object, X_blocks = X_blocks, weighted = TRUE)
        }
        if (consensus_block == "consensus") {
            block_object$variates[[consensus_block]] <-  .get_consensus_variates(object = block_object, X_blocks = X_blocks, weighted = FALSE)
        }
        block_object$names$blocks <- c(block_object$names$blocks, consensus_block)
        block_object$ncomp[consensus_block] <- block_object$ncomp[1]
        block_object$explained_variance[consensus_block] <- 0
        
    }
    
    block_object
}

## ------------------------------- mat.rank ------------------------------- ##
#' Matrix Rank
#' 
#' This function estimate the rank of a matrix.
#' 
#' \code{mat.rank} estimate the rank of a matrix by computing its singular
#' values \eqn{d[i]} (using \code{nipals}). The rank of the matrix can be
#' defined as the number of singular values \eqn{d[i] > 0}.
#' 
#' If \code{tol} is missing, it is given by
#' \code{tol=max(dim(mat))*max(d)*.Machine$double.eps}.
#' 
#' @param mat a numeric matrix or data frame that can contain missing values.
#' @param tol positive real, the tolerance for singular values, only those with
#' values larger than \code{tol} are considered non-zero.
#' @return The returned value is a list with components: \item{rank}{a integer
#' value, the matrix rank.} \item{tol}{the tolerance used for singular values.}
#' @author Sébastien Déjean, Ignacio González, Al J Abadi
#' @seealso \code{\link{nipals}}
#' @keywords algebra
#' @export
#' @examples
#' 
#' ## Hilbert matrix
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#' mat <- hilbert(16)
#' mat.rank(mat)
#' 
#' \dontrun{
#' ## Hilbert matrix with missing data
#' idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
#' m.na <- m <- hilbert(9)[, 1:6]
#' m.na[idx.na == 0] <- NA
#' mat.rank(m)
#' mat.rank(m.na)
#' }

mat.rank = function (mat, tol)
{
    
    if (length(dim(mat)) != 2) 
        stop("'mat' must be a numeric matrix.")
    
    mat = as.matrix(mat)
    
    if (!is.numeric(mat)) 
        stop("'mat' must be a numeric matrix.")
    
    d = nipals(mat)$eig
    max.d = d[1]
    min.d = d[length(d)]
    if (missing(tol)) 
        tol = max(dim(mat)) * max.d * .Machine$double.eps
    
    r = sum(d > tol)
    
    return(list(rank = r, tol = tol))
}

## ----------------------------- nearZeroVar ------------------------------ ##
#' Identification of zero- or near-zero variance predictors
#' 
#' Borrowed from the \pkg{caret} package. It is used as an internal function in
#' the PLS methods, but can also be used as an external function, in
#' particular when the data contain a lot of zeroes values and need to be
#' pre-filtered beforehand.
#' 
#' This function diagnoses predictors that have one unique value (i.e. are zero
#' variance predictors) or predictors that are have both of the following
#' characteristics: they have very few unique values relative to the number of
#' samples and the ratio of the frequency of the most common value to the
#' frequency of the second most common value is large.
#' 
#' For example, an example of near zero variance predictor is one that, for
#' 1000 samples, has two distinct values and 999 of them are a single value.
#' 
#' To be flagged, first the frequency of the most prevalent value over the
#' second most frequent value (called the ``frequency ratio'') must be above
#' \code{freqCut}. Secondly, the ``percent of unique values,'' the number of
#' unique values divided by the total number of samples (times 100), must also
#' be below \code{uniqueCut}.
#' 
#' In the above example, the frequency ratio is 999 and the unique value
#' percentage is 0.0001.
#' 
#' @param x a numeric vector or matrix, or a data frame with all numeric data.
#' @param freqCut the cutoff for the ratio of the most common value to the
#' second most common value.
#' @param uniqueCut the cutoff for the percentage of distinct values out of the
#' number of total samples.
#' @return \code{nearZeroVar} returns a list that contains the following
#' components:
#' 
#' \item{Position}{a vector of integers corresponding to the column positions
#' of the problematic predictors that will need to be removed.}
#' \item{Metrics}{a data frame containing the zero- or near-zero predictors
#' information with columns: \code{freqRatio}, the ratio of frequencies for the
#' most common value over the second most common value and,
#' \code{percentUnique}, the percentage of unique data points out of the total
#' number of data points.}
#' @author Max Kuhn, Allan Engelhardt, Florian Rohart, Benoit Gautier, AL J Abadi
#' for mixOmics
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}
#' @keywords utilities
#' @export
#' @examples
#' 
#' data(diverse.16S)
#' nzv = nearZeroVar(diverse.16S$data.raw)
#' length(nzv$Position) # those would be removed for the default frequency cut
nearZeroVar = function (x, freqCut = 95/5, uniqueCut = 10)
{
    
    if (is.vector(x))
        x = matrix(x, ncol = 1)
    
    freqRatio = apply(x, 2, function(data)
    {
        data = na.omit(data)
        
        if (length(unique(data)) == length(data))
        { # No duplicate
            return(1)
        } else if (length(unique(data)) == 1) { # Same value
            return(0)
        } else {
            t = table(data)
            return(max(t, na.rm = TRUE)/max(t[-which.max(t)], na.rm = TRUE))
        }
    })
    
    lunique = apply(x, 2, function(data) length(unique(data[!is.na(data)])))
    
    percentUnique = 100 * lunique/nrow(x)
    zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
    
    out = list()
    out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
    names(out$Position) = NULL
    out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
    out$Metrics = out$Metrics[out$Position, ]
    return(out)
}

## -------------------------- .check_test.keepX --------------------------- ##
#' check test.keepX
#'
#' @param test.keepX test.keepX
#' @param X X input from mixOmics tune models
#' @param indY indY
#' @param already.tested.X already.tested.X
#'
#' @return test.keepX, possibly re-ordered by names for list X
#' @noRd
#' @keywords Internal
#' @examples
#' 
.check_test.keepX <- function(test.keepX, 
                         X,
                         indY = NULL, # TODO
                         already.tested.X = NULL # TODO
                         )
{
    # TODO use this helper to check all test.keepX in the package
    ## -- checker for a pair of test.keepX and X
    .check_test.keepX_helper <- function(test.keepX_, X_)
    {
        if (is.data.frame(X_))
        {
            X_ <- as.matrix(X_)
        }
        
        if (mode(test.keepX_) != 'numeric' || !all(test.keepX_ %in% seq_len(ncol(X_))))
        {
            stop( "'test.keepX' values should be positive whole numbers < ncol(X) = ", 
                  ncol(X_), call. = FALSE)
        }
        invisible(NULL)
    }
    
    ## ------- X matrix ------- ##
    if (!is.null(dim(X)))
    {
        .check_test.keepX_helper(test.keepX, X)
    } 
    ## -------- X list -------- ##
    else
    {
        ## names
        if (! (is.list(test.keepX) && setequal(names(test.keepX), names(X)) ) )
        {
            stop("'test.keepX' must be a named list with names: names(X) = ", names(X))
        }
        else
        {
            test.keepX <- test.keepX[names(X)]
            mapply(z = test.keepX, w = X, FUN = function(z, w) .check_test.keepX_helper(z, w)) 
        }
        
    }
    invisible(test.keepX)
}
## ----------------------------- .check_ncomp ----------------------------- ##

## check that ncomp is positive integer <= smallest dim in the data
.check_ncomp <- function(ncomp, X)
{
    # TODO use this helper to check all ncomp in the package
    
    if (mode(ncomp) != 'numeric' || ncomp%%1 != 0)
        stop("'ncomp' must be a positive integer")
    ## ------- X matrix ------- ##
    if (!is.null(dim(X)))
    {
        X <- list(X)
    } 
    min.dim <- min(sapply(X, function(z) min(dim(z))))
    if (ncomp > min.dim)
        stop("'ncomp' must be smaller than or equal to the smallest dimenion in X: ", min.dim)
    
    invisible(ncomp)
}
