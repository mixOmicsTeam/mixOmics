## First example: analysis with sPLS-DA
data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)

# first tune on components only
tune = tune.splsda(X, Y, ncomp = 5, logratio = "none",
                   nrepeat = 10, folds = 10,
                   test.keepX = NULL, 
                   dist = "all",
                   progressBar = TRUE,
                   seed = 20) # set for reproducibility of example only
plot(tune) # optimal distance = centroids.dist
tune$choice.ncomp # optimal component number = 3

# then tune optimal keepX for each component
tune = tune.splsda(X, Y, ncomp = 3, logratio = "none",
                   nrepeat = 10, folds = 10, 
                   test.keepX = c(5, 10, 15), dist = "centroids.dist",
                   progressBar = TRUE,
                   seed = 20)

plot(tune)
tune$choice.keepX # optimal number of variables to keep c(15, 5, 15)

## With already tested variables:
tune = tune.splsda(X, Y, ncomp = 3, logratio = "none",
                   nrepeat = 10, folds = 10, 
                   test.keepX = c(5, 10, 15), already.tested.X = c(5, 10),
                   dist = "centroids.dist",
                   progressBar = TRUE,
                   seed = 20)
plot(tune)

## Second example: multilevel one-factor analysis with sPLS-DA

data(vac18)
X = vac18$genes
Y = vac18$stimulation
# sample indicates the repeated measurements
design = data.frame(sample = vac18$sample)

# tune on components
tune = tune.splsda(X, Y = Y, ncomp = 5, nrepeat = 10, logratio = "none",
                   test.keepX = NULL, folds = 10, dist = "max.dist", multilevel = design)

plot(tune)

# tune on variables
tune = tune.splsda(X, Y = Y, ncomp = 3, nrepeat = 10, logratio = "none",
                   test.keepX = c(5,50,100),folds = 10, dist = "max.dist", multilevel = design)

plot(tune)