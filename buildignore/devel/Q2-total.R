
data(nutrimouse)
X <- nutrimouse$gene  
Y <- nutrimouse$lipid

MyResult.pls <- pls(X,Y, ncomp = 4)  

set.seed(30)
## run it using RELEASE_3_10 & devel -- see why such a big discrepancy
perf.pls <- perf(MyResult.pls, validation = "Mfold", folds = 5,
                 progressBar = FALSE, nrepeat = 10)
perf.pls$measures$Q2.total

plot(perf.pls$Q2.total)

perf.pls$Q2.total

## these are not looking right in the vignette -- the optimal should have
## maximum cor

# tuning both X and Y
set.seed(30) # for reproducibility in this vignette, otherwise increase nrepeat
tune.spls.cor.XY <- tune.spls(X, Y, ncomp = 3,
                              test.keepX = c(8, 20, 50),
                              test.keepY = c(4, 8, 16),
                              validation = "Mfold", folds = 5,
                              nrepeat = 10, progressBar = FALSE,
                              measure = 'cor')
## visualise correlations
plot(tune.spls.cor.XY, measure = 'cor')
## visualise RSS
plot(tune.spls.cor.XY, measure = 'RSS')