## default method: shows cross correlation between 2 data sets
#------------------------------------------------------------------
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)

nutrimouse.sgccda <- block.splsda(X = data,
                                  Y = Y,
                                  design = "full",
                                  keepX = list(gene = c(10,10), lipid = c(15,15)),
                                  ncomp = 2)

cimDiablo(nutrimouse.sgccda, comp = c(1,2))
## change trim range
cimDiablo(nutrimouse.sgccda, comp = c(1,2), trim = 4)
## do not trim values
cimDiablo(nutrimouse.sgccda, comp = c(1,2), trim = FALSE)
