## Example: analysis with PLS-DA
data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)

# tune components and distance
tune = tune.plsda(X, Y, ncomp = 5, logratio = "none",
                   nrepeat = 10, folds = 10,
                   progressBar = TRUE,
                   seed = 20) # set for reproducibility of example only
plot(tune) # optimal distance = centroids.dist
tune$choice.ncomp # optimal component number = 3