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
.unexpected_err <- function(trying_to = NULL) {
    trying_to <- ifelse(is.null(trying_to), "", sprintf(" while trying to %s", trying_to))
    msg <- sprintf("Unexpected error%s. Please check the inputs and if problem persists submit an issue to https://github.com/mixOmicsTeam/mixOmics/issues", trying_to)
    stop(msg, call. = FALSE)
}
