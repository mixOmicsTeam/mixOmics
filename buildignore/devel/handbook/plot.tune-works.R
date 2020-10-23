
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

## Not run: 
set.seed(42)
tune.cor = tune.spls(X, Y, ncomp=2, test.keepX = c(5,10),test.keepY = c(5,8,10),, measure.tune = "cor", method = 'spls',
                 nrepeat=3, progressBar = TRUE, folds =3)
tune.cor
plot.tune.spls(tune.cor)

tune.RSS = tune.spls(X, Y, ncomp=2, test.keepX = c(5,10),test.keepY = c(5,8,10),, measure.tune = "RSS", method = 'spls',
                 nrepeat=3, progressBar = TRUE, folds =3)

# plot the results
plot.tune.spls(tune.RSS)
