## object of class 'spls'
# --------------------------
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50),
keepY = c(10, 10))

plotLoadings(toxicity.spls)

# with xlim
xlim = c(-0.5, 0.5)
plotLoadings(toxicity.spls, xlim = xlim)


\dontrun{
## object of class 'splsda'
# --------------------------
data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = as.factor(paste0('treatment_' ,liver.toxicity$treatment[, 4]))

splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# contribution on comp 1, based on the median.
# Colors indicate the group in which the median expression is maximal
plotLoadings(splsda.liver, comp = 1, method = 'median')
plotLoadings(splsda.liver, comp = 1, method = 'median', contrib = "max")

# contribution on comp 2, based on median.
#Colors indicate the group in which the median expression is maximal
plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = "max")

# contribution on comp 2, based on median.
# Colors indicate the group in which the median expression is minimal
plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = 'min')

# look at the contribution (median) for each variable
plot.contrib = plotLoadings(splsda.liver, comp = 2, method = 'median', plot = FALSE,
contrib = "max")
head(plot.contrib[[1]][,1:4])
# change the title of the legend and title name
plotLoadings(splsda.liver, comp = 2, method = 'median', legend.title = 'Time',
title = 'Contribution plot', contrib = "max")

# no legend
plotLoadings(splsda.liver, comp = 2, method = 'median', legend = FALSE, contrib = "max")

# change the color of the legend
plotLoadings(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4), contrib = "max")


# object 'splsda multilevel'
# -----------------

data(vac18)
X = vac18$genes
Y = vac18$stimulation
# sample indicates the repeated measurements
sample = vac18$sample
stimul = vac18$stimulation

# multilevel sPLS-DA model
res.1level = splsda(X, Y = stimul, ncomp = 3, multilevel = sample,
keepX = c(30, 137, 123))


name.var = vac18$tab.prob.gene[, 'Gene']
names(name.var) = colnames(X)

plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, size.name = 0.2, contrib = "max")

# too many transcripts? only output the top ones
plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, size.name = 0.5, ndisplay = 60, contrib = "max")



# object 'plsda'
# ----------------

# breast tumors
# ---
data(breast.tumors)
X = breast.tumors$gene.exp
Y = breast.tumors$sample$treatment

plsda.breast = plsda(X, Y, ncomp = 2)

name.var = as.character(breast.tumors$genes$name)
names(name.var) = colnames(X)

# with gene IDs, showing the top 60
plotLoadings(plsda.breast, contrib = 'max', comp = 1, method = 'median',
ndisplay = 60,
name.var = name.var,
size.name = 0.6,
legend.color = color.mixo(1:2))


# liver toxicity
# ---

data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$treatment[, 4]

plsda.liver = plsda(X, Y, ncomp = 2)

name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotLoadings(plsda.liver, contrib = 'max', comp = 1, method = 'median', ndisplay = 100,
name.var = name.var, size.name = 0.4,
legend.color = color.mixo(1:4))


# object 'sgccda'
# ----------------

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda = wrapper.sgccda(X = data,
Y = Y,
design = design,
keepX = list(gene = c(10,10), lipid = c(15,15)),
ncomp = 2)

plotLoadings(nutrimouse.sgccda,block=2)
plotLoadings(nutrimouse.sgccda,block="gene")

}


#' PCA example
#' -----------
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

#' Simple PCA plot
pca.liver = pca(X, ncomp = 2)
plotLoadings(pca.liver)

#' Customize PCA plot
plotLoadings(pca.liver, 
    comp = 2,
    ndisplay = 20,  # Show top 20 variables
    col = "steelblue",
    border = "black",
    size.name = 0.8,
    title = "PCA Loadings - Component 2",
    X.label = "Loading value",
    Y.label = "Variables")

#' PLS/SPLS example
#' ----------------
toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

#' Plot both blocks with custom layout
plotLoadings(toxicity.spls,
    comp = 2,
    block = c("X", "Y"),
    layout = c(2, 2),
    title = "PLS Loadings - Component 2",
    subtitle = c("Gene Block", "Clinical Block"))

#' PLSDA/SPLSDA example
#' -------------------
X = as.matrix(liver.toxicity$gene)
Y = as.factor(paste0('treatment_', liver.toxicity$treatment[, 4]))
splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))

#' Show contribution with gene names
name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    name.var = name.var,
    size.name = 0.5,
    title = "Liver Treatment - Component 2",
    legend.title = "Treatment",
    legend.color = color.mixo(1:4))


#' MINT PLSDA example
#' -----------------
data(stemcells)
X = stemcells$gene
Y = stemcells$celltype
study = stemcells$study
mint.splsda = mint.splsda(X = X, Y = Y, ncomp = 3, keepX = c(10, 5, 15), study = study)

#' All partial loadings with custom layout
plotLoadings(mint.splsda,
    comp = 2,
    study = "all.partial",
    contrib = "max",
    method = "median",
    layout = c(2, 4),
    legend.color = color.mixo(1:3),
    subtitle = paste("Study", 1:4),
    show.ties = TRUE,
    col.ties = "gray")

#' Advanced customization examples
#' -----------------------------

#' Custom variable names and sizes
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    name.var = name.var,
    size.name = 0.8,
    title = "Custom Variable Names",
    size.title = 1.2,
    size.axis = 0.8,
    size.labs = 1.0)

#' No legend
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    legend = FALSE)

#' Custom legend
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    legend.title = "Treatment Groups",
    legend.color = c("red", "blue", "green", "purple"),
    size.legend = 0.8)

#' Minimal display with top variables
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    ndisplay = 20,  # Show only top 20 variables
    title = "Top 20 Contributing Variables")

