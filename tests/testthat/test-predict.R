test_that("predict.mixo_pls works", code = {
    data("linnerud")
    X <- linnerud$exercise
    Y <- linnerud$physiological
    linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")
    
    indiv1 <- as.integer(colMeans(X)) + 3
    indiv2 <- as.integer(colMeans(X)) - 3
    newdata <- rbind(indiv1, indiv2)
    colnames(newdata) <- colnames(X)
    
    # yhat_approx <- colMeans(Y)
    
    pred <- predict(linn.pls, newdata)
    expect_equal(round(pred$predict[1]), round(175.8865 ))
    # pred$predict
    # yhat_approx
    
})


test_that("predict.mixo_plsda works", code = {
    data("liver.toxicity")
    X <- liver.toxicity$gene
    Y <- as.factor(liver.toxicity$treatment[, 4])
    
    
    ## if training is perfomed on 4/5th of the original data
    samp <- sample(1:5, nrow(X), replace = TRUE)
    test <- which(samp == 1)   # testing on the first fold
    train <- setdiff(1:nrow(X), test)
    
    plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
    test.predict <- predict(plsda.train, X[test, ], dist = "max.dist")
    expect_is(test.predict, "predict")
})

test_that("predict.block.splsda", code = {
    # example with block.splsda=diablo=sgccda and a missing block
    data(nutrimouse)
    # need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
    Y.mat = unmap(nutrimouse$diet)
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y.mat)
    # with this design, all blocks are connected
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
                    byrow = TRUE, dimnames = list(names(data), names(data)))
    
    # train on 75
    ind.train=NULL
    for(i in 1:nlevels(nutrimouse$diet))
        ind.train=c(ind.train,which(nutrimouse$diet==levels(nutrimouse$diet)[i])[1:6])
    
    #training set
    gene.train=nutrimouse$gene[ind.train,]
    lipid.train=nutrimouse$lipid[ind.train,]
    Y.mat.train=Y.mat[ind.train,]
    Y.train=nutrimouse$diet[ind.train]
    data.train=list(gene=gene.train,lipid=lipid.train,Y=Y.mat.train)
    
    #test set
    gene.test=nutrimouse$gene[-ind.train,]
    lipid.test=nutrimouse$lipid[-ind.train,]
    Y.mat.test=Y.mat[-ind.train,]
    Y.test=nutrimouse$diet[-ind.train]
    data.test=list(gene=gene.test,lipid=lipid.test)
    
    # example with block.splsda=diablo=sgccda and a missing block
    res.train = block.splsda(X=list(gene=gene.train,lipid=lipid.train),Y=Y.train,
                             ncomp=3,keepX=list(gene=c(10,10,10),lipid=c(5,5,5)))
    
    expect_warning((test.predict <- predict(res.train, newdata=data.test[2], method = "max.dist")))
    expect_is(test.predict$WeightedPredict, 'array')
})


test_that("predict.mint.splsda works", code = {
    ## example with mint.splsda
    data(stemcells)
    
    #training set
    ind.test = which(stemcells$study == "3")
    gene.train = stemcells$gene[-ind.test,]
    Y.train = stemcells$celltype[-ind.test]
    study.train = factor(stemcells$study[-ind.test])
    
    #test set
    gene.test = stemcells$gene[ind.test,]
    Y.test = stemcells$celltype[ind.test]
    study.test = factor(stemcells$study[ind.test])
    
    res = mint.splsda(X = gene.train, Y = Y.train, ncomp = 3, keepX = c(10, 5, 15),
                      study = study.train)
    
    pred = predict(res, newdata = gene.test, study.test = study.test)
    # saveRDS(pred, file = 'inst/tests/predict.mint.splsda.rda')
    pred_ref <- readRDS(system.file("testdata", "predict.mint.splsda.rda", package = 'mixOmics'))
    expect_identical(pred_ref, pred)
})
