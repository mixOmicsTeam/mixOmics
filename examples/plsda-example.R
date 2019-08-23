#' \dontrun{

data(breast.tumors)

## ---------------- with X as matrix and and Y as numeric/factor
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment
plsda.res1 <- plsda(X=X, Y=Y)
plotVar(plsda.res1, cutoff = 0.6)
plotLoadings(plsda.res1)

## ---------------- formula method for matrices
## 'formula' argument should be explicitly mentioned (formula = ...)
## for correct method dispatch
plsda.res2 <- plsda(formula = Y ~ X)
## exclude calls and see if all outputs  are identical
identical(plsda.res1[-1], plsda.res2[-1])
#> TRUE
## ---------------- MultiAssayExperiment and assay names as X and Y
## 'data' argument should be explicitly mentioned for correct method dispatch
plsda.res3 <- plsda(X='gene', Y='Dose.Group', data = liver.toxicity.mae)
plsda.res3

## ---------------- MultiAssayExperiment and formula with assay names
plsda.res4 <- plsda(formula = Dose.Group ~ gene, data = liver.toxicity.mae)
identical(plsda.res3[-1], plsda.res4[-1])
#> TRUE
#' }
