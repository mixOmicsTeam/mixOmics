library(mixOmicsData)
library(MultiAssayExperiment)
## This is a parametrised test set for combinations of X, Y, formula, and daya in pls
## the tests should satisfy the conditions explained in documentation of pls function.
## the 'miniACC' data with s
## ----------------------------------------------------- parameters to change
## MAE data
mae_data <- miniACC

## index of the X assay in MAE data
X_index <- 1 ## RNASeq2GeneNorm
## index of Y (when it's assay and not colData) assay in MAE data
Ya_index <- 2 ## gistict
## index of valid Y in colData(MAE data)
Yc_index <- 2 ## "years_to_birth" - must be numeric
## index of invalid numeric Y in colData(MAE data) (this will test informative error messages)
iYc_index_num <- 1 ## "patientID" - must be numeric
## index of invalid non-numeric Y in colData(MAE data) (this will test informative error messages)
iYc_index_char <- 7 ## pathologic_stage" - is in fact a valid class but not factor
## phenotypes column name for discriminatory analyses
phen_col <- "gender" ## uses Xa for single plsda and list(Xa, Ya) for block.plsda
## ----------------------------------------------------- no need to change the following

Xa <- names(assays(mae_data))[X_index] ## X assay
Ya <- names(assays(mae_data))[Ya_index] ## Y assay
f_Ya<- as.formula(paste0(Ya, " ~ ", Xa)) ## formula with Y assay

Yc <- names(colData(mae_data))[Yc_index] ## Y column data
f_Yc <- as.formula(paste0(Yc, " ~ ", Xa)) ## formula with Y column data

Xm_Yc <- t(as.matrix(assay(mae_data, Xa))) ## X matrix when Y column data
Xm_Ya <- t(as.matrix(assay(mae_data[,complete.cases(mae_data[,,c(Xa, Ya)])], Xa))) ## X matrix when Y is assay

Yam <- t(assay(mae_data[,complete.cases(mae_data[,,c(Xa, Ya)])], Ya)) ## Y assay matrix
Ycn <-  colData(mae_data[,,Xa])[,Yc] ## Y column data numeric

Y_inv <- colnames(colData(mae_data))[iYc_index_num] ## invalid coldata y
Y_inv.vec <- colData(mae_data)[,iYc_index_num] ## vector

Y_char <-colnames(colData(mae_data))[iYc_index_char]
Y_char.vec <- colData(mae_data)[,iYc_index_char]

phenXa <- colData(mae_data[,,Xa])[,phen_col] ## phenotype for plsda with Xa
phenXaYa <- colData(mae_data[,complete.cases(mae_data[,,c(Xa, Ya)])])[,phen_col] ## phenotype for block.plsda with Xa, Ya
