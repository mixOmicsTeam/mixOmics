# set up data
data(stemcells)
data <- stemcells$gene
type.id <- stemcells$celltype
exp <- stemcells$study

# tune number of components
tune_res <- tune.mint.plsda(X = data,Y = type.id, ncomp=5, 
                             near.zero.var=FALSE,
                             study=exp)

plot(tune_res)
tune_res$choice.ncomp # 1 component