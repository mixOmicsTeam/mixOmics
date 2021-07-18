data(liver.toxicity)

# implement IPCA on a microarray dataset
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")
ipca.res

# samples representation
plotIndiv(
    ipca.res,
    ind.names = as.character(liver.toxicity$treatment[, 4]),
    group = as.numeric(as.factor(liver.toxicity$treatment[, 4]))
)

\dontrun{
    plotIndiv(ipca.res,
              cex = 0.01,
              col = as.numeric(as.factor(liver.toxicity$treatment[, 4])),
              style = "3d")
}
# variables representation
plotVar(ipca.res, cex = 0.5)

\dontrun{
plotVar(ipca.res, rad.in = 0.5, cex = 0.5, style="3d")
}
