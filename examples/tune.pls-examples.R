# set up data
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

# tune PLS2 model to find optimal number of components
tune.res <- tune.pls(X, Y, ncomp = 10, measure = "cor",
                    folds = 5, nrepeat = 3, progressBar = TRUE)
plot(tune.res) # plot outputs

# PLS1 model example
Y1 <- liver.toxicity$clinic[,1]

tune.res <- tune.pls(X, Y1, ncomp = 10, measure = "cor",
                    folds = 5, nrepeat = 3, progressBar = TRUE)

plot(tune.res)

# Multilevel PLS2 model
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
design <- data.frame(sample = repeat.indiv)

tune.res <- tune.pls(X, Y1, ncomp = 10, measure = "cor", multilevel = design,
                     folds = 5, nrepeat = 3, progressBar = TRUE)

plot(tune)