#' Human Stem Cells Data
#' 
#' This data set contains the expression of a random subset of 400 genes in 125
#' samples from 4 independent studies and 3 cell types.
#' 
#' This data set contains the expression of a random subset of 400 genes in 125
#' samples from 4 independent studies and 3 cell types. Those studies can be
#' combined and analysed using the MINT procedure.
#' 
#' @name stemcells
#' @docType data
#' @usage data(stemcells)
#' @format A list containing the following components: \describe{
#' \item{list("gene")}{data matrix with 125 rows and 400 columns. Each row
#' represents an experimental sample, and each column a single gene.}
#' \item{list("celltype")}{a factor indicating the cell type of each sample.}
#' \item{list("study")}{a factor indicating the study from which the sample was
#' extracted.} }
#' @return none
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#' @keywords datasets
NULL
