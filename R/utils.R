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
        
}