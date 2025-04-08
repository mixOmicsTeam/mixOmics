
## 'spls' class - examples demonstrate how to control sample colours with sample names shown
# ------------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
spls.obj <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                 keepY = c(10, 10, 10))

# colors indicate time of necropsy, text is the dose, label at start of arrow
plotArrow(spls.obj,  group = as.factor(liver.toxicity$treatment[, 'Time.Group']),
          col = c("red", "blue", "purple", "darkgreen"), 
          ind.names  = liver.toxicity$treatment[, 'Dose.Group'],
          legend = TRUE, position.names = 'start', legend.title = 'Time.Group')


## 'rcc' class - examples demonstrate how to control shape of all samples
# -------------------------------------------------------------------------------

# create model
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
rcc.obj <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

primary_groups <- nutrimouse$diet

# plot samples coloured by primary groups, by default shapes are all circles
plotArrow(rcc.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE)

# plot samples coloured by primary groups, force all samples to have 
# a different shape (2 = triangle)
plotArrow(rcc.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE,
          pch = 2)


## 'sgccda' class - examples demonstrate how to control shape of different blocks
# --------------------------------------------------------------------------------

data(breast.TCGA)
idx = seq(1, length(breast.TCGA$data.train$subtype), 10)
X <- list(mRNA = breast.TCGA$data.train$mrna[idx,],
          miRNA = breast.TCGA$data.train$mirna[idx,],
          protein = breast.TCGA$data.train$protein[idx,])
Y <- breast.TCGA$data.train$subtype[idx] # set the response variable

diablo.obj <- block.splsda(X, Y, ncomp = 2) # undergo multiblock sPLS-DA

# plot the samples using an arrow plot - 
plotArrow(diablo.obj, 
          ind.names = FALSE,
          legend = TRUE,
          title = 'TCGA, DIABLO comp 1 - 2') 

pchs <- c(3, 2, 1)
names(pchs) <- c("miRNA", "mRNA", "protein")

# plot the samples using an arrow plot changed shapes of data blocks
plotArrow(diablo.obj, 
          ind.names = FALSE,
          legend = TRUE,
          title = 'TCGA, DIABLO comp 1 - 2',
          pch = pchs) 
