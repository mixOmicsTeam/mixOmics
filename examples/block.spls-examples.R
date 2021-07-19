# Example with multi omics TCGA study
# -----------------------------
data("breast.TCGA")
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna)
# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = c(2)
# set number of variables to select, per component and per data set (this is set arbitrarily)
list.keepX = list(mrna = rep(10, 2), mirna = rep(10,2))
list.keepY = c(rep(10, 2))

TCGA.block.spls = block.spls(X = data, Y = breast.TCGA$data.train$protein,
ncomp = ncomp, keepX = list.keepX, keepY = list.keepY, design = design)
TCGA.block.spls
# in plotindiv we color the samples per breast subtype group but the method is unsupervised!
plotIndiv(TCGA.block.spls, group =  breast.TCGA$data.train$subtype, ind.names = FALSE, legend=TRUE)
# illustrates coefficient weights in each block
plotLoadings(TCGA.block.spls, ncomp = 1)
plotVar(TCGA.block.spls, style = 'graphics', legend = TRUE)

## plot markers (selected markers) for mrna and mirna
group <- breast.TCGA$data.train$subtype
# mrna: show each selected feature separately and group by subtype
plotMarkers(object = TCGA.block.spls, comp = 1, block = 'mrna', group = group)
# mrna: aggregate all selected features, separate by loadings signs and group by subtype
plotMarkers(object = TCGA.block.spls, comp = 1, block = 'mrna', group = group, global = TRUE)
# proteins
plotMarkers(object = TCGA.block.spls, comp = 1, block = 'Y', group = group)
## only show boxplots
plotMarkers(object = TCGA.block.spls, comp = 1, block = 'Y', group = group, violin = FALSE)

\dontrun{
network(TCGA.block.spls)
}
