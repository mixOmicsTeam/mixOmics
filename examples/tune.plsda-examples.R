## Example: analysis with PLS-DA
data(breast.tumors)

# tune components and distance
tune = tune.plsda(breast.tumors$gene.exp, as.factor(breast.tumors$sample$treatment), 
                  ncomp = 5, logratio = "none",
                  nrepeat = 10, folds = 10,
                  progressBar = TRUE,
                  seed = 20) # set for reproducibility of example only
plot(tune) # optimal distance = centroids.dist
tune$choice.ncomp # optimal component number = 3

## Example: multilevel PLS-DA
data(vac18)
design <- data.frame(sample = vac18$sample) # set the multilevel design

tune1 <- tune.plsda(vac18$genes, vac18$stimulation, 
                    ncomp = 5, multilevel = design,
                   nrepeat = 10, folds = 10,
                   seed = 20)
plot(tune1)