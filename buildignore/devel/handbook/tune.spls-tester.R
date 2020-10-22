load('buildignore/devel/handbook/env.RData')
X <- X[,1:300]
Y <- Y[,1:600]
## ---- eval = FALSE-----------------------------------------------------------------------
list.keepX = c(seq(5, 50, 5), seq(60, 150, 10))
list.keepX = c(2, 5)
list.keepY = c(3:10)
list.keepY = c(3, 8)
ncomp = 2
nrepeat = 3


# update with tuning code later + outputs
set.seed(33)  # for reproducibility with this handbook, remove otherwise
sPLS.tune.reg.cor.liver = tune.spls(X, Y, folds= 3, test.keepX = list.keepX, test.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'regression', measure.tune = 'cor', method = 'spls')
PLS.tune.reg.cor.liver = tune.spls(X, Y, folds= 3, ncomp = ncomp, nrepeat = nrepeat, mode = 'regression', measure.tune = NULL, method = 'pls')
sPLS.tune.reg.cor.liver = tune.spls(X, Y, test.keepX = list.keepX, test.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'regression', measure.tune = 'cor', method = 'spls', progressBar = TRUE)

st <- list()
st$serial <- system.time({
    tune.spls(X, Y, folds= 30, ncomp = ncomp, nrepeat = 20, mode = 'regression', measure.tune = NULL, method = 'pls', BPPARAM = BiocParallel::SerialParam())
})

st$parallel <- system.time({
    tune.spls(X, Y, folds= 30, ncomp = ncomp, nrepeat = 20, mode = 'regression', measure.tune = NULL, method = 'pls', BPPARAM = BiocParallel::MulticoreParam(workers = 4))
})

st

# sPLS.tune.reg.cor.liver.ref <- sPLS.tune.reg.cor.liver
# PLS.tune.reg.cor.liver.ref <- PLS.tune.reg.cor.liver
# save(sPLS.tune.reg.cor.liver.ref, PLS.tune.reg.cor.liver.ref, file='res.RData')
load('res.RData')
# testthat::expect_equal(sPLS.tune.reg.cor.liver, sPLS.tune.reg.cor.liver.ref)
# testthat::expect_equal(PLS.tune.reg.cor.liver, PLS.tune.reg.cor.liver.ref)
waldo::compare(y = sPLS.tune.reg.cor.liver, sPLS.tune.reg.cor.liver.ref)
waldo::compare(y = PLS.tune.reg.cor.liver, PLS.tune.reg.cor.liver.ref)
