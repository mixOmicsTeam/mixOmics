context("tune.block.splsda")

test_that("tune.block.splsda works with and without parallel without auc", {

    data("breast.TCGA")
    data = list(
        mrna = breast.TCGA$data.train$mrna,
        mirna = breast.TCGA$data.train$mirna,
        protein = breast.TCGA$data.train$protein
    )
    design = matrix(
        1,
        ncol = length(data),
        nrow = length(data),
        dimnames = list(names(data), names(data))
    )
    diag(design) =  0
    ncomp <- 2
    nrep <- 3
    test.keepX = list(
        mrna = c(10, 20),
        mirna = c(20, 30),
        protein = c(3, 6)
    )
    
    ## -------------------- 1 cpu 
    set.seed(42)
    tune11 = tune.block.splsda(
        X = data,
        Y = breast.TCGA$data.train$subtype,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        nrepeat = nrep
    )
    expect_is(tune11, "tune.block.splsda")
    expect_equal(tune11$choice.ncomp$ncomp, 2)
    
    ## -------------------- parallel
    set.seed(42)
    tune41 = tune.block.splsda(
        X = data,
        Y = breast.TCGA$data.train$subtype,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        nrepeat = nrep,
        cpus = 2
    )
    expect_is(tune41, "tune.block.splsda")
    expect_equal(tune41$choice.ncomp$ncomp, 2)
    
    ## -------------------- already.tested.keepX
    already.tested.X = lapply(tune11$choice.keepX, function(x) {
        x[1]
    })
    
    tune42 = tune.block.splsda(
        X = data,
        Y = breast.TCGA$data.train$subtype,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        already.tested.X = already.tested.X,
        cpus = 2
    )
    
    expect_is(tune42, "tune.block.splsda")
    
})
