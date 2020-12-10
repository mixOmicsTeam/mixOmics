#' @section folds: 
#' During a cross-validation (CV), data are randomly split into \code{M}
#' subgroups (folds). \code{M-1} subgroups are then used to train submodels
#' which would be used to predict prediction accuracy statistics for the
#' held-out (test) data. All subgroups are used as the test data exactly once.
#' If \code{validation = "loo"}, leave-one-out CV is used where each group
#' consists of exactly one sample and hence \code{M == N} where N is the number
#' of samples.
