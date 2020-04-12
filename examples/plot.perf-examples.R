require(lattice)

## validation for objects of class 'pls' or 'spls'

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.pls <- pls(X, Y, ncomp = 3)
liver.perf <- perf(liver.pls, validation = "Mfold")

plot(liver.perf, criterion = "R2", layout = c(2, 2))

\dontrun{
## validation for objects of class 'plsda' or 'splsda'

data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
breast.perf <- perf(res, nrepeat = 5)

plot(breast.perf)
plot(breast.perf, col=1:3)
plot(breast.perf, col=1:3, sd=FALSE)

}
