
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

## Not run: 
set.seed(42)
tune.cor = tune.spls(X, Y, ncomp=2, test.keepX = c(2,4,8,12),test.keepY = c(5,6,7,8,9), measure.tune = "cor", method = 'spls',
                 nrepeat=3, progressBar = TRUE, folds =3, BPPARAM = BiocParallel::MulticoreParam(workers = 4, RNGseed = 23))
tune.cor
plot.tune.spls(tune.cor)
plot.tune.spls(tune.cor, pch.size = 3, cex = 1.8)
plot.tune.spls(tune.cor, measure = 'RSS')
plot.tune.spls(tune.cor, interactive = TRUE)

tune.RSS = tune.spls(X, Y, ncomp=2, test.keepX = c(5,10),test.keepY = c(5,8,10),, measure.tune = "RSS", method = 'spls',
                 nrepeat=3, progressBar = TRUE, folds =3, BPPARAM = BiocParallel::MulticoreParam(workers = 4))

# plot the results
plot.tune.spls(tune.RSS)
