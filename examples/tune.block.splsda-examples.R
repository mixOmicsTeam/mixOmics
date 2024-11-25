## Set up data

# load data
data("breast.TCGA")

# X data - list of mRNA and miRNA
X <- list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
          protein = breast.TCGA$data.train$protein)

# Y data - single data set of proteins
Y <- breast.TCGA$data.train$subtype

# subset the X and Y data to speed up computation in this example
set.seed(100)
subset <- mixOmics:::stratified.subsampling(breast.TCGA$data.train$subtype, folds = 3)[[1]][[1]]
X <- lapply(X, function(omic) omic[subset,])
Y <- Y[subset]

# set up a full design where every block is connected
# could also consider other weights, see our mixOmics manuscript
design = matrix(1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) =  0
design

## Tune number of components to keep
tune_res <- tune.block.splsda(X, Y, design = design,
                              ncomp = 5,
                              test.keepX = NULL,
                              nrepeat = 3,
                              seed = 13)

plot(tune_res)

tune_res$choice.ncomp # 3 components best

## Tune number of variables to keep

# definition of the keepX value to be tested for each block mRNA miRNA and protein
# names of test.keepX must match the names of 'data'
test.keepX = list(mrna = c(10, 30), mirna = c(15, 25), protein = c(4, 8))

# load parallel package
library(BiocParallel)

# run tuning in parallel on 2 cores
tune_res <- tune.block.splsda(X, Y, design = design,
                              ncomp = 2,
                              test.keepX = test.keepX,
                              nrepeat = 3,
                              seed = 13,
                              BPPARAM = SnowParam(workers = 2))

plot(tune_res)
tune_res$choice.keepX

# Now tuning a new component given previous tuned keepX
already.tested.X <- tune_res$choice.keepX
tune_res <- tune.block.splsda(X, Y, design = design,
                         ncomp = 3, 
                         test.keepX = test.keepX,
                         already.tested.X = already.tested.X,
                         BPPARAM = SnowParam(workers = 2))
tune_res$choice.keepX