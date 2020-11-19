# example with missing values where NIPALS is applied
# --------------------------------
data(multidrug)
X <- multidrug$ABC.trans
pca.res <- pca(X, ncomp = 4, scale = TRUE)
plot(pca.res)
print(pca.res)
biplot(pca.res, group = multidrug$cell.line$Class, legend.title = 'Class')

# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)))

# variable representation
plotVar(pca.res, var.names = TRUE, cutoff = 0.4, pch = 16)

\dontrun{
plotIndiv(pca.res, cex = 0.2,
    col = as.numeric(as.factor(multidrug$cell.line$Class)),style="3d")


plotVar(pca.res, rad.in = 0.5, cex = 0.5, style="3d")
}

# example with imputing the missing values using impute.nipals()
# --------------------------------
data("nutrimouse")
X <- data.matrix(nutrimouse$lipid)
X <- scale(X, center = TRUE, scale = TRUE)
## add missing values to X to impute and compare to actual values
set.seed(42)
na.ind <- sample(seq_along(X), size = 20)
true.values <- X[na.ind]
X[na.ind] <- NA
pca.no.impute <- pca(X, ncomp = 2)
plotIndiv(pca.no.impute, group = nutrimouse$diet, pch = 16)
X.impute <- impute.nipals(X, ncomp = 10)
## compare
cbind('imputed' = round(X.impute[na.ind], 2), 
      'actual' = round(true.values, 2))
## run pca using imputed matrix
pca.impute <- pca(X.impute, ncomp = 2)
plotIndiv(pca.impute, group = nutrimouse$diet, pch = 16)
# example with multilevel decomposition and CLR log ratio transformation 
# (ILR takes longer to run)
# ----------------
data("diverse.16S")
pca.res = pca(X = diverse.16S$data.TSS, ncomp = 3,
    logratio = 'CLR', multilevel = diverse.16S$sample)
plot(pca.res)
plotIndiv(pca.res, ind.names = FALSE, 
          group = diverse.16S$bodysite, 
          title = '16S diverse data',
          legend = TRUE, 
          legend.title = 'Bodysite')
