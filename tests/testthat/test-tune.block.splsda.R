context("tune.block.splsda")
library(BiocParallel)

test_that("tune.block.splsda works with and without parallel without auc", {

    data("breast.TCGA")
    data = list(
        mrna = breast.TCGA$data.train$mrna,
        mirna = breast.TCGA$data.train$mirna,
        protein = breast.TCGA$data.train$protein
    )
    
    design = 'full'
    ncomp <- 2
    folds <- 3
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
    tune11 = tune.block.splsda(
        X = data,
        Y = Y,
        folds = folds,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        nrepeat = nrep,
        BPPARAM = SerialParam(), seed = 42
    )
    expect_is(tune11, "tune.block.splsda")
    
    ## -------------------- parallel
    tune41 = tune.block.splsda(
        X = data,
        Y = Y,
        folds = folds,
        ncomp = ncomp,
        test.keepX = test.keepX,
        design = design,
        nrepeat = nrep,
        BPPARAM = SnowParam(workers = 2), seed = 42
    )
    expect_equal(tune11$choice.keepX, tune41$choice.keepX)
    
    # ## -------------------- already.tested.keepX
    # already.tested.X = lapply(tune11$choice.keepX, function(x) {
    #     x[1]
    # })
    # 
    # set.seed(42)
    # tune42 = tune.block.splsda(
    #     X = data,
    #     Y = Y,
    #     ncomp = ncomp,
    #     folds = folds,
    #     test.keepX = test.keepX,
    #     design = design,
    #     already.tested.X = already.tested.X,
    #     BPPARAM = SerialParam(RNGseed = 42)
    # )
    # expect_equal(tune11$choice.keepX, tune42$choice.keepX)
    
})

test_that("(tune.block.splsda:error): catches invalid values of 'folds'", {
    
    data("breast.TCGA")
    
    samples <- c(1:3, 50:52, 79:81)
    data = list(miRNA = breast.TCGA$data.train$mirna[samples,], 
                mRNA = breast.TCGA$data.train$mrna[samples,], 
                proteomics = breast.TCGA$data.train$protein[samples,])
    Y = breast.TCGA$data.train$subtype[samples]

    design = matrix(0.1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
    diag(design) = 0 # set diagonal to 0s

    # set grid of values for each component to test
    test.keepX = list (mRNA = c(1,2), 
                       miRNA = c(1,2), 
                       proteomics = c(1,2))
    
    expect_error(tune.block.splsda(X = data, Y = Y, ncomp = 2, 
                              test.keepX = test.keepX, design = design, folds=10),
                 "'folds' cannot be greater than the number of input samples",
                 fixed=T)
    
    expect_error(tune.block.splsda(X = data, Y = Y, ncomp = 2, 
                              test.keepX = test.keepX, design = design, folds=1),
                 "'folds' needs to be at least 2",
                 fixed=T)
    
    expect_error(tune.block.splsda(X = data, Y = Y, ncomp = 2, 
                              test.keepX = test.keepX, design = design, folds="random.value"),
                 "'folds' need to be non-NULL and numeric",
                 fixed=T)
})

test_that("tune.block.splsda works independently and in tune wrapper the same", {
  
  data("breast.TCGA")
  data = list(
    mrna = breast.TCGA$data.train$mrna,
    mirna = breast.TCGA$data.train$mirna,
    protein = breast.TCGA$data.train$protein
  )
  
  design = 'full'
  ncomp <- 2
  folds <- 3
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
    folds = folds,
    ncomp = ncomp,
    test.keepX = test.keepX,
    design = design,
    nrepeat = nrep,
    BPPARAM = SerialParam(), seed = 42
  )
  expect_is(tune11, "tune.block.splsda")
  
  ## -------------------- parallel
  set.seed(42)
  tune41 = tune(
    method = "block.splsda",
    X = data,
    Y = Y,
    folds = folds,
    ncomp = ncomp,
    test.keepX = test.keepX,
    design = design,
    nrepeat = nrep,
    BPPARAM = SerialParam(), seed = 42
  )
  expect_equal(tune11$choice.keepX, tune41$choice.keepX)
  
  # check can plot outputs
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(tune41))
  expect_silent(plot(tune11))

})

test_that("tune.block.splsda works when test.keepX = NULL and gives same result as perf()", code = {
  
  data("breast.TCGA")
  data = list(
    mrna = breast.TCGA$data.train$mrna,
    mirna = breast.TCGA$data.train$mirna,
    protein = breast.TCGA$data.train$protein
  )
  
  design = 'full'
  ncomp <- 2
  folds <- 3
  nrep <- 3
  test.keepX = list(
    mrna = c(10, 20),
    mirna = c(20, 30),
    protein = c(3, 6)
  )
  
  set.seed(100)
  subset <- mixOmics:::stratified.subsampling(breast.TCGA$data.train$subtype, folds = 4)[[1]][[1]]
  X <- lapply(data, function(omic) omic[subset,])
  Y <- breast.TCGA$data.train$subtype[subset]
  
  # tune on components only
  tune_res <- suppressWarnings(
    tune.block.splsda(X, Y, ncomp = 3, design = design,
              nrepeat = 1, folds = 3,
              test.keepX = NULL, seed = 20)
  )
  
  # run perf
  block.splsda_res <- block.splsda(X, Y, ncomp = 3, design = design)
  perf_res <- suppressWarnings(
    perf(block.splsda_res, ncomp = 3, nrepeat = 1, folds = 3, dist = "max.dist", seed = 20)
  )
  
  # check results match
  expect_equal(tune_res$error.rate$mrna[1,1], perf_res$error.rate$mrna[1,1])
})
