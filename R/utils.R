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
## ----------- .check_cpus ----------- 
#' Check cpus argument
#'
#' @param cpus Integer, the number of cpus for parallel processing.
#'
#' @return The number of cpus available/applicable, or condition.
#'
#' @noRd
.check_cpus <- function(cpus) {
    if (!is.numeric(cpus) || cpus <= 0)
        stop("Number of CPUs should be a positive integer.", call. = FALSE)
    
    if (cpus < 0 || is.null(cpus))
        cpus <- 1
    
    if (cpus > detectCores())
        message(sprintf("\nOnly %s CPUs available for parallel processing.\n", detectCores()))
    return(cpus)
        
}

## ----------- .check_alpha ----------- 
#' Check significance threshold sanity
#'
#' @param alpha numeric, significance threshold for t.test
#'
#' @return NULL if acceptable value, otherwise condition
#'
#' @noRd
.check_alpha <- function(alpha=NULL) {
    if (is.null(alpha))
        alpha <- 0.01
    
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1)
        stop("invalid 'signif.threshold'. Use 0.01 or 0.05.", call. = FALSE)
    alpha
}

## ----------- .unexpected_err ----------- 
#' Unexpected error handler for the package
#'
#' To be used in unexperimented situations where handlers fail
#' @param trying_to character, the context in which unexpected error occurred.
#'
#' @return Error message
#'
#' @noRd
.unexpected_err <- function(trying_to = NULL) {
    trying_to <- ifelse(is.null(trying_to), "", sprintf(" while trying to %s", trying_to))
    msg <- sprintf("Unexpected error%s. Please check the inputs and if problem persists submit an issue to https://github.com/mixOmicsTeam/mixOmics/issues", trying_to)
    stop(msg, call. = FALSE)
}

## ----------- .on_unix ----------- 
#' Check OS type for parallel processing
#'
#' @return Logical, FALSE if windows OS, TRUE if unix OS.
#'
#' @noRd
.onUnix <- function() {
    return(ifelse(.Platform$OS.type == "unix", TRUE, FALSE))
}

## ----------- .unlist_repeat_cv_output ----------- 
#' repeat_cv_perf.diablo helper to unlist internal outputs to
#' previous list format for downstream work
#'
#' @return list of outputs needed for downstream work
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

## ----------- stratified_subsampling ----------- 
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

## ----------- .name_list ----------- 
#' Create a named list
#'
#' Creates a named list of a character vector where names and values
#' are the same, useful for lapply family.
#' @param char Character vector.
#'
#' @return A named list.
#'
#' @noRd
.name_list <- function(char) {
    out <- as.list(char)
    names(out) <- char
    return(out)
}

## ----------- .mixo_rng ----------- 
#' RNG used for package tests
#'
#' Creating a function so that it can easily be changed.
#' It's mainly used in form of the given example in unit tests.
#'
#' @return An R version.
#' @example RNGversion(.mixo_rng())
#' @noRd
.mixo_rng <- function() {"3.6.0"}

## ----------- .change_if_null ----------- 
#' Set default value if NULL
#'
#' For arguments with default NULL value, sets the desired default value
#' 
#' @param arg The function argument with NULL/no default value
#' @param default The desired default value
#'
#' @return The desired default value
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

## ----------- .add_consensus_blocks ----------- 
#' Add consensus blocks to DIABLO object
#'
#' For plotIndiv(blocks = "consensus")
#' 
#' @param block_object A diablo object
#' @param consensus_blocks One or both of c('consensus', 'weighted.consensus')
#'
#' @return A diablo object with consensus blocks added for smooth input into plotIndiv.sgccda
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
