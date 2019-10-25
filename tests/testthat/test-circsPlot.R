test_that("circosPlot works", code = {
    data(nutrimouse)
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    
    nutrimouse.sgccda <- wrapper.sgccda(X=data,
                                        Y = Y,
                                        design = design,
                                        keepX = list(gene=c(8,8), lipid=c(4,4)),
                                        ncomp = 2,
                                        scheme = "horst")
    
    
    cp_res <- circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1,
                        color.Y = 1:5, color.blocks = c("green","brown"), color.cor = c("magenta", "purple"))
    expect_is(cp_res, "matrix")
})