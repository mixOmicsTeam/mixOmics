pca.assay <- 'lipid'
# formals(sipca)[c('ncomp', 'mode', 'keepX')] <- list(3,'deflation', c(50,50,50))
datasets <- data(package="mixOmicsData")$results[,"Item"]
data(list=datasets)
