

## Arrow plot of individuals for objects of class 'spls'
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


## individual representation for objects of class 'sgccda'
# ----------------------------------------------------

data(nutrimouse)
Y <- nutrimouse$diet
data <- list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 <- matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)
sgccda.obj <- 
    wrapper.sgccda(X = data,
                   Y = Y,
                   design = design,
                   ncomp = 2,
                   keepX = list(gene = c(10,10), lipid = c(15,15)))


# by default colours correspond to outcome Y, shapes to block
plotArrow(sgccda.obj)

# change pch - TO DO





