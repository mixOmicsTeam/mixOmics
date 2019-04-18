context("test-auroc")

test_that("auroc works", {
    data(breast.tumors)
    set.seed(1)
    test=sample(1:47,5,replace=FALSE)
    X <- breast.tumors$gene.exp
    X.test<-breast.tumors$gene.exp[test,]
    Y <- breast.tumors$sample$treatment
    Y.test<-breast.tumors$sample$treatment[test]
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    auc.plsda=auroc(res.plsda,plot = TRUE,roc.comp = 1)
    
    expect_equal(matrix(auc.plsda$Comp1),
    rbind(0.863, 2.473e-05))

    expect_equal(matrix(auc.plsda$Comp2),
    rbind(0.9981, 7.124e-09))

 })
