context("plotLoadings")

test_that("plotLoadings.spls works", code = {
    data(liver.toxicity)
    X = liver.toxicity$gene
    Y = liver.toxicity$clinic
    
    toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50),
                         keepY = c(10, 10))
    
    pl_res <- plotLoadings(toxicity.spls)
    
    expect_is(pl_res, "list")
    
})

test_that("plotLoadings.splsda works", code = {
    data(liver.toxicity)
    X = as.matrix(liver.toxicity$gene)
    Y = as.factor(liver.toxicity$treatment[, 4])
    
    splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))

    pl_res <- plotLoadings(splsda.liver, comp = 1, method = 'median')
    
    expect_is(pl_res, "list")
    
})

test_that("plotLoadings.block.splsda works", code = {
    data(nutrimouse)
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutrimouse.sgccda = wrapper.sgccda(X = data,
                                       Y = Y,
                                       design = design,
                                       keepX = list(gene = c(10,10), lipid = c(15,15)),
                                       ncomp = 2)
    
    pl_res <- plotLoadings(nutrimouse.sgccda,block=2)
    
    expect_is(pl_res, "list")
    
})

test_that("plotLoadings.mint.splsda works", code = {
    data(stemcells)
    data = stemcells$gene
    type.id = stemcells$celltype
    exp = stemcells$study
    
    res = mint.splsda(X = data, Y = type.id, ncomp = 3, keepX = c(10,5,15), study = exp)
    pl_res <- plotLoadings(res, contrib = "max")
    
    expect_is(pl_res, "list")
    
})


test_that("plotLoadings.mint.spls works", code = {
    data(stemcells)
    samples <- c(1:5,60:64)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$gene[samples+5, 1:10]
    S <- as.character(stemcells$study[samples])
    
    res = mint.spls(X = X, Y = Y, ncomp = 3, 
                    keepX = seq(3, 9, 3), 
                    keepY = seq(3, 9, 3), 
                    study = S)
    pl_res <- plotLoadings(res, contrib = "max")
    
    expect_is(pl_res, "list")
})


test_that("plotLoadings margin errrors is handled properly", code = {
    data(nutrimouse)
    Y = nutrimouse$diet
    gene = nutrimouse$gene
    lipid = nutrimouse$lipid
    ## extend feature names
    suff <- "-a-long-suffix-from-abolutely-nowhere-which-is-gonna-be-longer-than-margins-and-you-best-believe-that-this-suffix-is-too-long"
    colnames(gene) <- paste0(colnames(gene), suff)
    colnames(lipid) <- paste0(colnames(lipid), suff)
    data = list(gene = gene, lipid = lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutrimouse.sgccda = block.splsda(X = data,
                                     Y = Y,
                                     design = design,
                                     keepX = list(gene = c(10,10), lipid = c(15,15)),
                                     ncomp = 2)
    expect_error(plotLoadings(nutrimouse.sgccda, contrib = "min"), regexp = "plotLoadings encountered margin errors")
})
