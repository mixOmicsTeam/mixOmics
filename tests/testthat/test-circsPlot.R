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

test_that("circosPlot works with similar feature names in different blocks", code = {

    create_similar_feature_names <- function(data_list)
    {
        lapply(data_list, function(x){
            colnames(x) <- paste0('feature_', seq_len(ncol(x)))
            x
        })
    }
    
    
    
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, 
                mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    
    data <- create_similar_feature_names(data)
    list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = rep(10, 2))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                     ncomp = 2, keepX = list.keepX, design = 'full')
    cp_res <- .quiet(circosPlot(TCGA.block.splsda, cutoff = 0.7))
    expect_is(cp_res, "matrix")
})

test_that("circosPlot works when using the indY parameter", code = {
    
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, 
                mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    
    list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = rep(10, 2))
    TCGA.block.spls = block.spls(X = data, indY = 3, 
                                     ncomp = 2, keepX = list.keepX, design = 'full')
    cp_res <- circosPlot(TCGA.block.spls, cutoff = 0.7, group = breast.TCGA$data.train$subtype)
    
    expect_is(cp_res, "matrix")
})

