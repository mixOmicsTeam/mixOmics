# block.splsda
# -------------
data("breast.TCGA")
# this is the X data as a list of mRNA, miRNA and proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
protein = breast.TCGA$data.train$protein)
# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = c(2)
# set number of variables to select, per component and per data set (this is set arbitrarily)
list.keepX = list(mrna = rep(5,2), mirna = rep(5,2), protein = rep(5,2))


TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                 ncomp = ncomp, keepX = list.keepX, design = design)
## use design = 'full'
TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                 ncomp = ncomp, keepX = list.keepX, design = 'full')
TCGA.block.splsda$design

plotIndiv(TCGA.block.splsda, ind.names = FALSE)
## use design = 'null'
TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                 ncomp = ncomp, keepX = list.keepX, design = 'null')
TCGA.block.splsda$design
## set all off-diagonal elements to 0.5
TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                 ncomp = ncomp, keepX = list.keepX, design = 0.5)
TCGA.block.splsda$design
# illustrates coefficient weights in each block
plotLoadings(TCGA.block.splsda, ncomp = 1, contrib = 'max')
plotVar(TCGA.block.splsda, style = 'graphics', legend = TRUE)

## plot markers (selected variables) for mrna and mirna
# mrna: show each selected feature separately
plotMarkers(object = TCGA.block.splsda, comp = 1, block = 'mrna')
# mrna: aggregate all selected features and separate by loadings signs
plotMarkers(object = TCGA.block.splsda, comp = 1, block = 'mrna', global = TRUE)
# proteins
plotMarkers(object = TCGA.block.splsda, comp = 1, block = 'protein')
# show top 5 markers
plotMarkers(object = TCGA.block.splsda, comp = 1, block = 'protein', markers = 1:5)
# show specific markers
my.markers <- selectVar(TCGA.block.splsda, comp = 1)[['protein']]$name[c(1,3,5)]
my.markers
plotMarkers(object = TCGA.block.splsda, comp = 1, block = 'protein', markers = my.markers)
