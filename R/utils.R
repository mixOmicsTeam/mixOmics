## ----------- .trim_long_names ----------- 
#' Trim long variable names
#'
#' Trims feature names longer than \code{len} adding \code{...} at the end of
#' the name with a message. This helps avoid margin errors in plots. 
#' \code{len=23} is chosen as default only because longest ensemble id for
#' mouse/human genes/transcripts is this length.
#' 
#' @param var_names character vector of variable names
#' @param len integer, names longer than this will be trimmed so total length
#' will be \code{len}, including \code{...}.
#'
#' @return character vector of trimmed names
#'
#' @examples
#' .trim_long_names(var_names = c("long-variable-name-of-length-31", "short-variable"))
#' @noRd
.trim_long_names <- function(var_names=NULL, len=23) {
    if (any(nchar(var_names) > len)) {
        message(sprintf("Some variable names are too long. Trimmed for visualisation purposes."))
                
        var_names <- sapply(var_names, function(char) {
            ifelse(nchar(char) > len, sprintf("%s...", substring(char, 1, len-3)) , char)
        })
    }
    
    return(var_names)
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
