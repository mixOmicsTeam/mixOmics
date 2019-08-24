library(mixOmicsData)
library(MultiAssayExperiment)

## This is a parametrised test set for combinations of X, Y, formula, and daya in pls
## the tests should satisfy the conditions explained in documentation of pls function.

## -------------- parameters to change --------------
## MAE data
block.mae <- breast.TCGA.mae

## index of the X assays in MAE data
X_indices <- 1:2 ## mrna, mirna
## index of Y (when it's assay and not colData) assay in MAE data
Y_index <- 3 ## protein
## index of valid Y in colData(MAE data)
Yc_index <- 1 ## "subtype" - must be numeric

## -------------- no need to change the following --------------
## the names used down here must be unique across all helper-* files

X_bpls <- as.list(experiments(block.mae)[X_indices])
Y_bpls <- experiments(block.mae)[[Y_index]]

## transpose
X_bpls <- lapply(X_bpls, t)
Y_bpls <- t(Y_bpls)
