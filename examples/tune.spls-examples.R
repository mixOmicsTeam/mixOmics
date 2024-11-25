# set up data
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

# tune spls model for components only
tune.res.ncomp <- tune.spls( X, Y, ncomp = 5,
                             test.keepX = NULL,
                             test.keepY = NULL, measure = "cor",
                             folds = 5, nrepeat = 3, progressBar = TRUE)
plot(tune.res.ncomp) # plot outputs

# tune spls model for number of X and Y variables to keep
tune.res <- tune.spls( X, Y, ncomp = 3,
                      test.keepX = c(5, 10, 15),
                      test.keepY = c(3, 6, 8), measure = "cor",
                      folds = 5, nrepeat = 3, progressBar = TRUE)
plot(tune.res) # plot outputs