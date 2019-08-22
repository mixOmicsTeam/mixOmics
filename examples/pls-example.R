#' \dontrun{
library(mixOmicsData)
data(linnerud)
## ---------------- with X and Y as matrices
X <- linnerud$exercise
Y <- linnerud$physiological
pls.res1 <- pls(X=X, Y=Y)
plotVar(pls.res1)

## ---------------- formula method for matrices
## 'formula' argument should be explicitly mentioned (formula = ...)
## for correct method dispatch
pls.res2 <- pls(formula = Y ~ X)
## exclude calls and see if all outputs  are identical
identical(pls.res1[-1], pls.res2[-1])
#> TRUE
## ---------------- MultiAssayExperiment and assay names as X and Y
## 'data' argument should be explicitly mentioned for correct method dispatch
data("linnerud.mae")
pls.res3 <- pls(X='exercise', Y='physiological', data = linnerud.mae, mode = "regression")
identical(pls.res1[-1], pls.res3[-1])
#> TRUE

## ---------------- MultiAssayExperiment and formula with assay names
pls.res4 <- pls(formula = physiological ~ exercise, data = linnerud.mae)
identical(pls.res1[-1], pls.res4[-1])
#> TRUE

## ---------------- data=MultiAssayExperiment; X=assay and Y=colData
data("liver.toxicity")
data("liver.toxicity.mae")
toxicity.pls1 <- pls(data = liver.toxicity.mae,  formula = Dose.Group~gene, ncomp = 3)
toxicity.pls2 <- pls(data = liver.toxicity.mae,  Y='Dose.Group', X='gene', ncomp = 3)
identical(toxicity.pls1[-1], toxicity.pls2[-1])
#> TRUE

#' }
