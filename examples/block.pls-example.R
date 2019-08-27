#' \dontrun{
# Example with TCGA multi omics study
# -----------------------------------
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
library(mixOmicsData)
data("breast.TCGA")
protein <- breast.TCGA$data.train$protein
mrna <- breast.TCGA$data.train$mrna
mirna <- breast.TCGA$data.train$mirna

# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = c(2)

TCGA.block.pls = block.pls(X = list(mrna=mrna, mirna=mirna), Y = protein, ncomp = ncomp,
design = design)

## alternatively, we can use a formula notation
TCGA.block.pls.formula = block.pls(formula = protein ~ mrna + mirna, ncomp = ncomp,
                           design = design)
## check if the results will be identical
identical(TCGA.block.pls, TCGA.block.pls.formula)
#> TRUE
# in plotindiv we color the samples per breast subtype group but the method is unsupervised!
# here Y is the protein data set
plotIndiv(TCGA.block.pls, group =  breast.TCGA$data.train$subtype, ind.names = FALSE)

# Example with MultiAssayExperiments
# -----------------------------------
library(mixOmicsData)
data("breast.TCGA.mae")

TCGA.block.pls.mae <- block.pls(data = breast.TCGA.mae, formula = protein ~ mrna + mirna, ncomp = ncomp,
                             design = design)
## check if the results will be identical
identical(TCGA.block.pls, TCGA.block.pls.mae)
#> TRUE

#' }