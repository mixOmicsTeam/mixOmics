#' @section nrepeat: 
#' The cross-validation process is repeated \code{nrepeat} times and the
#' accuracy measures are averaged across repeats. If \code{validation = "loo"},
#' the process does not need to be repeated as there is only one way to split N
#' samples into N groups and hence nrepeat is forced to be 1.
