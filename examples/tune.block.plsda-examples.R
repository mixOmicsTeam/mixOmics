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
tune_res <- tune.block.plsda(X, Y, design = design,
                              ncomp = 5,
                              nrepeat = 3,
                              seed = 13,
                             dist = c("max.dist", "centroids.dist"))

plot(tune_res)

tune_res$choice.ncomp # 3 components best for max.dist, 1 for centroids.dist


## Tune number of components to keep - use weighted vote rather than majority vote
tune_res <- tune.block.plsda(X, Y, design = design,
                             ncomp = 5,
                             nrepeat = 3,
                             seed = 13,
                             dist = c("max.dist", "centroids.dist"),
                             weighted = FALSE)
plot(tune_res)

tune_res$weights