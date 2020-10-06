
## plot of individuals for objects with two datasets only (X and Y)
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

## plot of individuals for objects of class 'pls' or 'spls'
# ----------------------------------------------------
plotArrow(nutri.res)
## customise the ggplot object as you wish
plotArrow(nutri.res) + geom_vline(xintercept = 0, alpha = 0.5) + 
    geom_hline(yintercept = 0, alpha = 0.5) +
    labs(x = 'Dim 1' , y = 'Dim 2', title = 'Nutrimouse') +
    theme_minimal() +
    xlim(-4, 4)
## individual name position
plotArrow(nutri.res, ind.names.position = 'end')
plotArrow(nutri.res, scale.blocks = FALSE) ## do not scale blocks to same X scales
plotArrow(nutri.res, comp = c(1,3))
## custom pch
plotArrow(nutri.res, pch = 10, pch.size = 3)
plotArrow(nutri.res, pch = c(X = 1, Y = 0))
## custom arrow
plotArrow(nutri.res, arrow.alpha = 0.6, arrow.size = 0.6, arrow.length = 0.15)

## group samples
plotArrow(nutri.res, group = nutrimouse$genotype)
plotArrow(nutri.res, group = nutrimouse$genotype, legend.title = 'Genotype')

## custom ind.names
plotArrow(nutri.res,
           ind.names = paste0('ID', rownames(nutrimouse$gene)), 
           ind.names.size = 3)

## plot of individuals for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                      keepY = c(10, 10, 10))

# colors indicate time of necropsy, text is the dose, label at start of arrow
plotArrow(toxicity.spls,  group = liver.toxicity$treatment[, 'Time.Group'],
           ind.names  = liver.toxicity$treatment[, 'Dose.Group'],
           legend = TRUE, position.names = 'start', legend.title = 'Time.Group')

## individual representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca <- wrapper.sgcca(X = data,
                                  design = design1,
                                  penalty = c(0.3, 0.5, 1),
                                  ncomp = 3,
                                  scheme = "centroid")

plotArrow(nutrimouse.sgcca, group = nutrimouse$genotype, ind.names = TRUE, 
           legend.title = 'Genotype' )

## custom pch by block
blocks <- names(nutrimouse.sgcca$variates)
pch <- seq_along(blocks)
names(pch) <- blocks
pch
#>   gene   lipid     Y 
#>   1       2        3 
plotArrow(nutrimouse.sgcca, group = nutrimouse$genotype, ind.names = TRUE, 
           pch =pch, legend.title = 'Genotype' )

## individual representation for objects of class 'sgccda'
# ----------------------------------------------------
# Note: the code differs from above as we use a 'supervised' GCCA analysis
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- 
    wrapper.sgccda(X = data,
                   Y = Y,
                   design = design1,
                   ncomp = 2,
                   keepX = list(gene = c(10,10), lipid = c(15,15)),
                   scheme = "centroid")


## Default colours correspond to outcome Y
plotArrow(nutrimouse.sgccda1)

