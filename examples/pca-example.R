#' \dontrun{

library(mixOmics.data)
# example with missing values where NIPALS is applied
# --------------------------------
pca.res <- pca(X = multidrug$ABC.trans, ncomp = 4, scale = TRUE)
pca.res
plot(pca.res)
print(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

## example with MultiAssayExperiment class
## --------------------------------
pca.res <- pca(data = multidrug.mae, X = 'ABC.trans', ncomp = 4, scale = TRUE)
plot(pca.res)
print(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

## samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
          group = as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, cex = 0.2,
            col = as.numeric(as.factor(multidrug$cell.line$Class)),style="3d")


## variable representation
plotVar(pca.res)
## 3D
plotVar(pca.res, rad.in = 0.5, cex = 0.5,style="3d")


## example with multilevel decomposition and CLR log ratio transformation
## (ILR longer to run)
## ----------------
pca.res = pca(X = diverse.16S$data.TSS, ncomp = 5,
              logratio = 'CLR', multilevel = diverse.16S$sample)
plot(pca.res)
plotIndiv(pca.res, ind.names = FALSE, group = diverse.16S$bodysite,
          title = '16S diverse data', legend = TRUE)

#' }
