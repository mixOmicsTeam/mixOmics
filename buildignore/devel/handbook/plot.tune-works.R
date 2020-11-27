
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

set.seed(42)
## -----  measure = 'cor'
tune.cor = tune.spls(X, Y, ncomp=2, test.keepX = c(2,4),test.keepY = c(5,10),measure.tune = 'cor',
                     nrepeat=3, progressBar = TRUE, folds =3)
## outputs
names(tune.cor)
# print
tune.cor
## plot
plot(tune.cor, measure = 'cor')
plot(tune.cor, measure = 'RSS')
tune.cor$choice.keepX
tune.cor$choice.keepY

## -----  measure = 'RSS'

tune.RSS = tune.spls(X, Y, ncomp=2, test.keepX = c(2,4),test.keepY = c(5,10),measure.tune = 'RSS',
                     nrepeat=3, progressBar = TRUE, folds =3)
## automatically detects measure
plot(tune.RSS)
tune.RSS$choice.keepX
tune.RSS$choice.keepY
