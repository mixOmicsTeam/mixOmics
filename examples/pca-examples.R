# example with missing values where NIPALS is applied
# --------------------------------
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plot(pca.res)
print(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

\dontrun{
# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
group = as.numeric(as.factor(multidrug$cell.line$Class)))

\dontrun{
plotIndiv(pca.res, cex = 0.2,
col = as.numeric(as.factor(multidrug$cell.line$Class)),style="3d")
}

# variable representation
plotVar(pca.res)

\dontrun{
plotVar(pca.res, rad.in = 0.5, cex = 0.5,style="3d")
}

# example with multilevel decomposition and CLR log ratio transformation (ILR longer to run)
# ----------------

data("diverse.16S")
pca.res = pca(X = diverse.16S$data.TSS, ncomp = 5,
logratio = 'CLR', multilevel = diverse.16S$sample)
plot(pca.res)
plotIndiv(pca.res, ind.names = FALSE, group = diverse.16S$bodysite, title = '16S diverse data',
legend = TRUE)
}