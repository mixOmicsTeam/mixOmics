# set up data
data(stemcells)
data <- stemcells$gene
type.id <- stemcells$celltype
exp <- stemcells$study

# tune number of components
tune_res <- tune.mint.splsda(X = data,Y = type.id, ncomp=5, 
                             near.zero.var=FALSE,
                             study=exp,
                             test.keepX = NULL)

plot(tune_res)
tune_res$choice.ncomp # 1 component

## tune number of variables to keep
tune_res <- tune.mint.splsda(X = data,Y = type.id, ncomp = 1, 
                             near.zero.var = FALSE,
                             study=exp,
                             test.keepX=seq(1,10,1))

plot(tune_res)
tune_res$choice.keepX # 9 variables to keep on component 1

## only tune component 3 and keeping 10 genes on comp1
tune_res <- tune.mint.splsda(X = data, Y = type.id, ncomp = 2, study = exp,
                             already.tested.X = c(9),
                             test.keepX = seq(1,10,1))
plot(tune_res)
tune_res$choice.keepX # 10 variables to keep on comp2