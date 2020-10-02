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
    
    set.seed(100)
    subset <- mixOmics:::stratified.subsampling(breast.TCGA$data.train$subtype, folds = 4)[[1]][[1]]
    data <- lapply(data, function(omic) omic[subset,])
    Y <- breast.TCGA$data.train$subtype[subset]
    
    ## -------------------- 1 cpu 
    set.seed(42)
    tune11 = tune.block.splsda(
        X = data,
        Y = Y,
        folds = 3,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        nrepeat = nrep
    )
    expect_is(tune11, "tune.block.splsda")
    expect_equal(tune11$choice.keepX, list(mrna = c(10, 10), mirna = c(20, 20), protein = c(3, 3)))
    
    ## -------------------- parallel
    BPPARAM <- if (!.onUnix()) BiocParallel::SnowParam(workers = 2) else BiocParallel::MulticoreParam(workers = 2)
    tune41 = tune.block.splsda(
        X = data,
        Y = Y,
        folds = 3,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        nrepeat = nrep,
        BPPARAM = BPPARAM
    )
    expect_is(tune41, "tune.block.splsda")
    
    ## -------------------- already.tested.keepX
    already.tested.X = lapply(tune11$choice.keepX, function(x) {
        x[1]
    })
    
    tune42 = tune.block.splsda(
        X = data,
        Y = Y,
        ncomp = ncomp,
        folds = 3,
        test.keepX = test.keepX,
        design = design,
        already.tested.X = already.tested.X,
        BPPARAM = BPPARAM
    )
    expect_is(tune42, "tune.block.splsda")
    
})
