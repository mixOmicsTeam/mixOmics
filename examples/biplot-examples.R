data("nutrimouse")
## --------- pca ---------- ##
pca.lipid <- pca(nutrimouse$lipid, ncomp = 3, scale = TRUE)
# seed for reproducible geom_text_repel

data("nutrimouse")
## --------- pca ---------- ##
pca.lipid <- pca(nutrimouse$lipid, ncomp = 3, scale = TRUE)
# seed for reproducible geom_text_repel
set.seed(42)
biplot(pca.lipid, group = as.factor(c(rep("A", 20), rep("B", 20))))

biplot(pca.lipid, group = as.factor(c(rep("A", 20), rep("B", 20))),
       col = c("red", "blue"))

set.seed(42)
biplot(pca.lipid)
## correlation cutoff to filter features
biplot(pca.lipid, cutoff = c(0.8))
## tailor threshold for each component
biplot(pca.lipid, cutoff = c(0.8, 0.7))
## customise components
biplot(pca.lipid, cutoff = c(0.8), comp = c(1,3))

## customise ggplot in an arbitrary way
biplot(pca.lipid) + theme_linedraw() + 
    # add vline
    geom_vline(xintercept = 0, col = 'green') +
    # add hline
    geom_hline(yintercept = 0, col = 'green') +
    # customise labs
    labs(x = 'Principal Component 1', y = 'Principal Component 2')

## group samples
biplot(pca.lipid, group = nutrimouse$diet, legend.title = 'Diet')

## customise variable labels
biplot(pca.lipid, 
       var.names.col = color.mixo(2),
       var.names.size = 4,
       var.names.angle = TRUE
)

## no arrows
biplot(pca.lipid, group = nutrimouse$diet, legend.title = 'Diet', 
       var.arrow.col = NULL, var.names.col = 'black')

## add x=0 and y=0 lines in function
biplot(pca.lipid, group = nutrimouse$diet, legend.title = 'Diet', 
       var.arrow.col = NULL, var.names.col = 'black', 
       vline = TRUE, hline = TRUE)

## --------- spca
## example with spca
spca.lipid <- spca(nutrimouse$lipid, ncomp = 2, scale = TRUE, keepX = c(8, 6))
biplot(spca.lipid, var.names.col = 'black', group = nutrimouse$diet, 
       legend.title = 'Diet')

## --------- pls ---------- ##
data("nutrimouse")
pls.nutrimouse <- pls(X = nutrimouse$gene, Y = nutrimouse$lipid, ncomp = 2)
biplot(pls.nutrimouse, group = nutrimouse$genotype, block = 'X',
       legend.title = 'Genotype', cutoff = 0.878) 

biplot(pls.nutrimouse, group = nutrimouse$genotype, block = 'Y',
       legend.title = 'Genotype', cutoff = 0.8) 

## --------- plsda ---------- ##
data(breast.tumors)
X <- breast.tumors$gene.exp
colnames(X) <- paste0('GENE_', colnames(X))
rownames(X) <- paste0('SAMPLE_', rownames(X))
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
biplot(plsda.breast, cutoff = 0.72)
## remove arrows
biplot(plsda.breast, cutoff = 0.72, var.arrow.col = NULL, var.names.size = 4)
