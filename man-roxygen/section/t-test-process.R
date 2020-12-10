#' @section t-test-process: 
#' The optimisation process is data-driven and similar to the process detailed
#' in (Rohart et al., 2016), where one-sided t-tests assess whether there is a
#' gain in performance when incrementing the number of features or components in
#' the model. However, it will assess all the provided grid through pair-wise
#' comparisons as the performance criteria do not always change linearly with
#' respect to the added number of features or components.
