#' \dontrun{

data(liver.toxicity)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
print(spca.rat)
## variable representation
plotVar(spca.rat, cex=2)
## 3d plot
plotVar(spca.rat,style="3d")
# }

# example with MultiAssayExperiment class
# --------------------------------
## X is the name of an assay from data
spca.rat <- spca(data=liver.toxicity.mae, X='gene', ncomp = 3, keepX = rep(50, 3))
plotVar(spca.rat, cex=2)
# \dontrun{


## samples representation
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3],
          group = as.numeric(liver.toxicity$treatment[, 3]))
plotIndiv(spca.rat, cex = 0.01,
          col = as.numeric(liver.toxicity$treatment[, 3]),style="3d")


# example with multilevel decomposition and CLR log ratio transformation
# ----------------

data("diverse.16S")
spca.res = spca(X = diverse.16S$data.TSS, ncomp = 5, keepX = rep(20, 5),
              logratio = 'CLR', multilevel = diverse.16S$sample)
plotIndiv(spca.res, ind.names = FALSE, group = diverse.16S$bodysite, title = '16S diverse data - sPCA',
          legend=TRUE)

  #' }
