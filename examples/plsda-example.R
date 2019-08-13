#' \dontrun{
library(mixOmics.data)

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
plsda.res3 <- plsda(X='exercise', Y='physiological', data = linnerud.mae)
identical(plsda.res1[-1], plsda.res3[-1])
#> TRUE

## ---------------- MultiAssayExperiment and formula with assay names
plsda.res4 <- plsda(formula = physiological ~ exercise, data = linnerud.mae, mode = "classic",)
identical(plsda.res1[-1], plsda.res4[-1])
#> TRUE

## ---------------- MultiAssayExperiment; X=assay and Y=colData
toxicity.plsda1 <- plsda(data = liver.toxicity.mae,  formula = Dose.Group~gene, ncomp = 3)
toxicity.plsda2 <- plsda(data = liver.toxicity.mae,  Y='Dose.Group', X='gene', ncomp = 3)
identical(toxicity.plsda1[-1], toxicity.plsda2[-1])
#> TRUE

#' }
