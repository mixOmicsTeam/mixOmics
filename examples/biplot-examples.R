data("mtcars")
pca.mtcars <- pca(mtcars, ncomp = 2, scale = TRUE)
# seed for reproducible geom_text_repel
set.seed(42)
biplot(pca.mtcars)

## cutomise ggplot in an arbitrary way
biplot(pca.mtcars) + theme_linedraw() + 
    # add vline
    geom_vline(xintercept = 0, col = 'green') +
    # add hline
    geom_hline(yintercept = 0, col = 'green') +
    # customise labs
    labs(x = 'Principal Component 1', y = 'Principal Component 2')

## group samples
biplot(pca.mtcars, group = mtcars$cyl, legend.title = 'Cyl')

## customise variable labels
biplot(pca.mtcars, 
       var.names.col = color.mixo(2),
       var.names.size = 4,
       var.names.angle = TRUE
       )

## no arrows
biplot(pca.mtcars, group = mtcars$cyl, legend.title = 'Cyl', 
       var.arrow.col = NULL, var.names.col = 'black')

## add x=0 and y=0 lines in function
biplot(pca.mtcars, group = mtcars$cyl, legend.title = 'Cyl', 
       var.arrow.col = NULL, var.names.col = 'black', 
       vline = TRUE, hline = TRUE)

## example with spca
spca.mtcars <- spca(mtcars, ncomp = 2, scale = TRUE, keepX = c(8, 6))
biplot(spca.mtcars, var.names.col = 'black', group = mtcars$gear, 
       legend.title = 'Gear')
