\dontrun{
data("breast.TCGA")
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
            protein = breast.TCGA$data.train$protein)
# set up a full design where every block is connected
# could also consider other weights, see our mixOmics manuscript
design = matrix(1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = 3

# Tuning the first two components
# -------------
## Not run: 
# definition of the keepX value to be tested for each block mRNA miRNA and protein
# names of test.keepX must match the names of 'data'
test.keepX = list(mrna = c(10, 30), mirna = c(15, 25), protein = c(4, 8))

# the following may take some time to run, so we subset the data first.
# Note that for thorough tuning, nrepeat should be >= 3 so that significance of 
# the model improvement can be measured
## ---- subset by 3rd of samples
set.seed(100)
subset <- mixOmics:::stratified.subsampling(breast.TCGA$data.train$subtype, folds = 3)[[1]][[1]]
data <- lapply(data, function(omic) omic[subset,])
Y <- breast.TCGA$data.train$subtype[subset]
## ---- run
tune <- tune.block.splsda(
    X = data,
    Y = Y,
    ncomp = ncomp,
    test.keepX = test.keepX,
    design = design,
    nrepeat = 2, 
    BPPARAM = MulticoreParam(workers = detectCores()-1)
)

plot(tune)
tune$choice.ncomp
tune$choice.keepX

# Now tuning a new component given previous tuned keepX
already.tested.X = tune$choice.keepX
tune = tune.block.splsda(X = data, Y = Y,
                         ncomp = 4, test.keepX = test.keepX, design = design,
                         already.tested.X = already.tested.X,
                         BPPARAM = MulticoreParam(workers = detectCores()-1)
                         )
tune$choice.keepX
}