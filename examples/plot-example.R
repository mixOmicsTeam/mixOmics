#' \dontrun{
## successful: TRUE
## final: TRUE

library(mixOmics.data)
pca.res <- pca(multidrug$ABC.trans, ncomp = 10, scale = TRUE)
plot(pca.res)
#' }
