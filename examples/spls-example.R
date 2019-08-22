#' \dontrun{
## 1st example: accepted input formats
#--------------------------------------------------------------
data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
keepX <- c(50, 50)
keepY <- c(5, 5)
spls.res1 <- spls(X=X, Y=Y, keepX = keepX, keepY = keepY)
plotVar(spls.res1)

## ---------------- formula method for matrices
## 'formula' argument should be explicitly mentioned (formula = ...)
## for correct method dispatch
spls.res2 <- spls(formula = Y ~ X, keepX = keepX, keepY = keepY)
## exclude calls and see if all outputs  are identical
identical(spls.res1[-1], spls.res2[-1])
#> TRUE
## ---------------- MultiAssayExperiment and assay names as X and Y
data(liver.toxicity.mae)
## 'data' argument should be explicitly mentioned for correct method dispatch
spls.res3 <- spls(X='gene', Y='clinic', data = liver.toxicity.mae,
                  keepX = keepX, keepY = keepY)
identical(spls.res1[-1], spls.res3[-1])
#> TRUE

## ---------------- MultiAssayExperiment and formula with assay names
spls.res4 <- spls(formula = clinic ~ gene, data = liver.toxicity.mae,
                  keepX = keepX, keepY = keepY)
identical(spls.res1[-1], spls.res4[-1])
#> TRUE

## ---------------- MultiAssayExperiment; X=assay and Y=colData
toxicity.spls1 <- spls(data = liver.toxicity.mae, formula = Dose.Group~gene,
                       ncomp = 2, keepX=keepX)
toxicity.spls2 <- spls(data = liver.toxicity.mae, Y='Dose.Group', X='gene',
                       ncomp = 2, keepX=keepX)
identical(toxicity.spls1[-1], toxicity.spls2[-1])
#> TRUE

## 2nd example: one-factor multilevel analysis with sPLS, selecting a subset of variables
#--------------------------------------------------------------

# note: we made up those data, pretending they are repeated measurements
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each

# this is a spls (unsupervised analysis) so no need to mention any factor in design
# we only perform a one level variation split
design <- data.frame(sample = repeat.indiv)
res.spls.1level <- spls(X = liver.toxicity$gene,
                        Y=liver.toxicity$clinic,
                        multilevel = design,
                        ncomp = 3,
                        keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                        mode = 'canonical')

# set up colors and pch for plotIndiv
col.stimu <- 1:nlevels(design$stimu)

plotIndiv(res.spls.1level, rep.space = 'X-variate', ind.names = FALSE,
          group = liver.toxicity$treatment$Dose.Group,
          pch = 20, main = 'Gene expression subspace',
          legend = TRUE)


plotIndiv(res.spls.1level, rep.space = 'Y-variate', ind.names = FALSE,
          group = liver.toxicity$treatment$Dose.Group,
          pch = 20, main = 'Clinical measurements ssubpace',
          legend = TRUE)

plotIndiv(res.spls.1level, rep.space = 'XY-variate', ind.names = FALSE,
          group = liver.toxicity$treatment$Dose.Group,
          pch = 20, main = 'Both Gene expression and Clinical subspaces',
          legend = TRUE)

## 3rd example: two-factor multilevel analysis with sPLS, selecting a subset of variables
#--------------------------------------------------------------

dose <- as.factor(liver.toxicity$treatment$Dose.Group)
time <- as.factor(liver.toxicity$treatment$Time.Group)
# note: we made up those data, pretending they are repeated measurements
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
design <- data.frame(sample = repeat.indiv, dose = dose, time = time)

res.spls.2level = spls(liver.toxicity$gene,
                       Y = liver.toxicity$clinic,
                       multilevel = design,
                       ncomp=2,
                       keepX = c(10,10), keepY = c(5,5))

#' }
