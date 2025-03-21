## 'pca' class - examples demonstrate how to control sample colours and shapes
# -----------------------------------------------------------------------------

# subset data and create model
data("srbct")
X <- srbct$gene[1:6, ]
rownames(X) <- c(paste0("Sample_", 1:6))
pca.obj <- pca(X, ncomp = 3)

primary_groups <- as.factor(c(rep("Group_1", 2), rep("Group_2", 2), rep("Group_3", 2)))
# [1] Group_1 Group_1 Group_2 Group_2 Group_3 Group_3
# Levels: Group_1 Group_2 Group_3
secondary_groups <- as.factor(c(rep("A", 3), rep("B", 2), rep("C", 1)))
# [1] A A A B B C
# Levels: A B C

# plot samples coloured by primary groups, show sample names
plotIndiv(pca.obj, ind.names = TRUE,
          group = primary_groups, legend = TRUE)

# plot samples coloured using custom colours by primary groups, show sample names
plotIndiv(pca.obj, ind.names = TRUE,
          group = primary_groups, legend = TRUE,
          col = c("red", "pink", "blue"))

# plot samples coloured by primary groups, by default shapes match primary groups
plotIndiv(pca.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE)

# plot samples coloured by primary groups, force all samples to have the same shape (2 = triangle)
plotIndiv(pca.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE,
          pch = 2)

# plot samples coloured by primary groups, use shapes to visualise secondary grouping
plotIndiv(pca.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE,
          pch = secondary_groups)

# plot samples coloured by primary groups, use shapes to visualise secondary grouping 
# and change order of secondary groups
plotIndiv(pca.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE,
          pch = factor(secondary_groups, levels = c("B", "C", "A")))


## 'rcc' class - examples demonstrate how to control rep.space
# ------------------------------------------------------------

# create model
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
rcc.obj <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# plot samples, by default makes a panel plot for X and Y subspaces (multi)
plotIndiv(rcc.obj)

# plot samples only on X-variate subspace
plotIndiv(rcc.obj, rep.space = "X-variate")

# plot samples only on XY-variate subspace
plotIndiv(rcc.obj, rep.space = "XY-variate")


## 'spls' class - examples demonstrate how to add ellipses and centroids/stars on groups
# --------------------------------------------------------------------------------------

# create model
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
spls.obj <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                      keepY = c(10, 10, 10))

# plot samples with ellipse on groups
plotIndiv(spls.obj, group = liver.toxicity$treatment$Time.Group, ellipse = TRUE)

# plot samples with centroids on groups
plotIndiv(spls.obj, group = liver.toxicity$treatment$Time.Group, centroid = TRUE)

# plot samples with centroids and stars on groups
plotIndiv(spls.obj, group = liver.toxicity$treatment$Time.Group, centroid = TRUE, star = TRUE)



## 'splsda' class - examples demonstrate how to add ellipses and backgrounds based on 
# predicted classes
# -----------------------------------------------------------------------------------------------------

# create model
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment
splsda.obj <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# plot samples with ellipse on groups, note groups do not have to be defined
plotIndiv(splsda.obj, ellipse = TRUE, ellipse.level = 0.8)

# plot samples with background coloured by predicted classes
background <- background.predict(splsda.obj, comp.predicted = 2, dist = "max.dist")
plotIndiv(splsda.obj, background = background)



## 'sgccda' class - examples demonstrate how to control which data blocks are plotted
# ------------------------------------------------------------------------------------

# create model
data(nutrimouse)
Y <- nutrimouse$diet
data <- list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design <- matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)
sgccda.obj <- wrapper.sgccda(X = data, Y = Y, design = design, ncomp = 2,
                             keepX = list(gene = c(10,10), lipid = c(15,15)))

# plot samples, by default one data block for each plot
plotIndiv(sgccda.obj)

# plot samples for just the gene data block
plotIndiv(sgccda.obj, blocks = 1)
plotIndiv(sgccda.obj, blocks = "gene")

# plot samples by averaging components from all blocks
plotIndiv(sgccda.obj, blocks = "average")

# plot samples by the weighted average of the components according to their correlation with Y
plotIndiv(sgccda.obj, blocks = "weighted.average")



## 'mint.splsda' class - examples demonstrate how to control which studies are plotted
# ------------------------------------------------------------------------------------

# create model
data(stemcells)
mint.obj <- mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2, 
                  keepX = c(10, 5), study = stemcells$study)

# plot samples, by default samples are plotted together coloured by groups and pch by study
plotIndiv(mint.obj, legend = TRUE)

# plot samples separated by study, can control layout
plotIndiv(mint.obj, legend = TRUE, study = "all.partial")
plotIndiv(mint.obj, legend = TRUE, study = "all.partial", layout = c(1,1))