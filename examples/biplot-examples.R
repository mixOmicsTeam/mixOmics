data("nutrimouse")
pca.lipid <- pca(nutrimouse$lipid, ncomp = 2, scale = TRUE)
# seed for reproducible geom_text_repel
set.seed(42)
biplot(pca.lipid)
## correlation cutoff to filter features
biplot(pca.lipid, cutoff = c(0.8))
## tailor threshold for each component
biplot(pca.lipid, cutoff = c(0.8, 0.7))

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

## example with spca
spca.lipid <- spca(nutrimouse$lipid, ncomp = 2, scale = TRUE, keepX = c(8, 6))
biplot(spca.lipid, var.names.col = 'black', group = nutrimouse$diet, 
       legend.title = 'Diet')
