#' \dontrun{
## successful: TRUE
## final: TRUE

data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 10, scale = TRUE)
plot(pca.res)
#' }
