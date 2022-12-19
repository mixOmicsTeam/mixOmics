context("auroc")

test_that("auroc works", {
    
    data(breast.tumors)
    set.seed(1)
    test=sample(1:47,5,replace=FALSE)
    X <- breast.tumors$gene.exp
    X.test<-breast.tumors$gene.exp[test,]
    Y <- breast.tumors$sample$treatment
    Y.test<-breast.tumors$sample$treatment[test]
    
    res.plsda <- plsda(X, Y, ncomp = 2)
    
    auc.plsda=.quiet(auroc(res.plsda,plot = TRUE,roc.comp = 1))
    
    expect_equal(matrix(auc.plsda$Comp1),
    rbind(0.863, 2.473e-05))

    expect_equal(matrix(auc.plsda$Comp2),
    rbind(0.9981, 7.124e-09))

 })


test_that("Safely handles zero var (non-zero center) features", {
  
  X1 <- data.frame(matrix(rnorm(100000, 5, 5), nrow = 100))
  X2 <- data.frame(matrix(rnorm(150000, 5, 5), nrow = 100))
  Y <- c(rep("A", 50), rep("B", 50))
  
  X <- list(block1=X1, block2=X2)
  
  list.keepX <- list(block1=c(15, 15), block2=c(30,30))
  
  set.seed(9425)
  X$block1[,1] <- rep(1, 100)
  model = suppressWarnings(block.splsda(X = X, Y = Y, ncomp = 2,
                                        keepX = list.keepX, design = "full",
                                        near.zero.var = T))
  
  auc.splsda = .quiet(auroc(model))
  
  .expect_numerically_close(auc.splsda$block1$comp1[[1]], 0.815)
  .expect_numerically_close(auc.splsda$block2$comp2[[2]], 2.22e-16)
  
})


test_that("combined auroc", {
  
  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna,
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)
  Y <- breast.TCGA$data.train$subtype
  
  model1 <- splsda(X$miRNA, Y, ncomp=3)
  model2 <- splsda(X$mRNA, Y, ncomp=3)
  model3 <- splsda(X$proteomics, Y, ncomp=3)
  
  models <- list(miRNA.Model=model1,
            mRNA.Model=model2,
            proteomics.Model=model3) 
  
  auroc.out <- auroc(models, print = F, plot=F, roc.comp = 3)
  
  .expect_numerically_close(auroc.out$auc$miRNA.Model[1],
                            0.993)
  
  .expect_numerically_close(auroc.out$auc$proteomics.Model[1],
                            0.9945)
})

dev.off()
unlink(list.files(pattern = "*.pdf"))
