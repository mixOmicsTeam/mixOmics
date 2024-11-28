## sPLS2 model example (more than one Y outcome variable)

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

## sPLS1 model example (only one Y outcome variable)

# set up data 
Y1 <- liver.toxicity$clinic[,1]

# tune spls model for components only
plot(tune.spls(X, Y1, ncomp = 3, 
               folds = 3, 
               test.keepX = NULL, test.keepY = NULL))

# tune spls model for number of X variables to keep, note for sPLS1 models 'measure' needs to be set
plot(tune.spls(X, Y1, ncomp = 3, 
               folds = 3, measure = "MSE",
               test.keepX = c(5, 10, 15), test.keepY = c(3, 6, 8)))


## sPLS2 multilevel model example

# set up multilevel design
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
design <- data.frame(sample = repeat.indiv)

# tune spls model for components only
tune.res.ncomp <- tune.spls( X, Y, ncomp = 5,
                             test.keepX = NULL,
                             test.keepY = NULL, measure = "cor", multilevel = design,
                             folds = 5, nrepeat = 3, progressBar = TRUE)
plot(tune.res.ncomp) # plot outputs

# tune spls model for number of X and Y variables to keep
tune.res <- tune.spls( X, Y, ncomp = 3,
                       test.keepX = c(5, 10, 15),
                       test.keepY = c(3, 6, 8), measure = "cor", multilevel = design,
                       folds = 5, nrepeat = 3, progressBar = TRUE)
plot(tune.res) # plot outputs
