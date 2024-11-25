# set up data
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

# tune pls model to find optimal number of components
tune.res <- tune.pls( X, Y, ncomp = 10, measure = "cor",
                             folds = 5, nrepeat = 3, progressBar = TRUE)
plot(tune.res) # plot outputs

# PLS1 model
Y1 <- liver.toxicity$clinic[,1]

tune.res <- tune.pls( X, Y1, ncomp = 10, measure = "cor",
                      folds = 5, nrepeat = 3, progressBar = TRUE)

plot(tune.res)