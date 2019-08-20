#' \dontrun{

data(liver.toxicity)

# implement IPCA on a microarray dataset
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode = "deflation")

# samples representation
plotIndiv(
  ipca.res,
  ind.names = as.character(liver.toxicity$treatment[, 4]),
  group = as.numeric(as.factor(liver.toxicity$treatment[, 4]))
)

# example with MultiAssayExperiment class
# --------------------------------

ipca.res <-
  ipca(liver.toxicity.mae,
       assay = 'gene',
       ncomp = 3,
       mode = "deflation")
ipca.res


# variables representation with cutoff
plotVar(ipca.res, cex = 1, cutoff = 0.5)

## 3d
plotVar(
  ipca.res,
  rad.in = 0.5,
  cex = 0.5,
  style = "3d",
  cutoff = 0.8
)

  #' }
