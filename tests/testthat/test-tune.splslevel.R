context("tune.splslevel")
library(BiocParallel)

test_that("tune.splslevel works and is the same in parallel and when run in tune wrapper", code = {
  
  # set up data
  data(liver.toxicity)
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  design <- data.frame(sample = repeat.indiv)
  
  # run in serial
  set.seed(42)
  tune.splslevel.res.1<- tune.splslevel(X = liver.toxicity$gene,
                                      Y=liver.toxicity$clinic,
                                      multilevel = design,
                                      test.keepX = c(5,10,15),
                                      test.keepY = c(1,2,5),
                                      ncomp = 1,
                                      BPPARAM = SerialParam(RNGseed = 42))
  
  # run in parallel
  set.seed(42)
  tune.splslevel.res.2<- tune.splslevel(X = liver.toxicity$gene,
                                        Y=liver.toxicity$clinic,
                                        multilevel = design,
                                        test.keepX = c(5,10,15),
                                        test.keepY = c(1,2,5),
                                        ncomp = 1,
                                        BPPARAM = SnowParam(RNGseed = 42, workers = 2))
  
  # in tune wrapper in serial
  set.seed(42)
  tune.splslevel.res.3<- tune(X = liver.toxicity$gene,
                                        Y=liver.toxicity$clinic,
                                        multilevel = design,
                                        test.keepX = c(5,10,15),
                                        test.keepY = c(1,2,5),
                                        ncomp = 1,
                                        BPPARAM = SerialParam(RNGseed = 42),
                              method = "spls")
  
  # in tune wrapper in parallel
  set.seed(42)
  tune.splslevel.res.4<- tune(X = liver.toxicity$gene,
                              Y=liver.toxicity$clinic,
                              multilevel = design,
                              test.keepX = c(5,10,15),
                              test.keepY = c(1,2,5),
                              ncomp = 1,
                              BPPARAM = SnowParam(RNGseed = 42, workers = 2),
                              method = "spls")
  
  
  # check outputs
  .expect_numerically_close(tune.splslevel.res.1$cor.value[1,1], 0.9637933)
  .expect_numerically_close(tune.splslevel.res.2$cor.value[1,1], 0.9637933)
  .expect_numerically_close(tune.splslevel.res.3$cor.value[1,1], 0.9637933)
  .expect_numerically_close(tune.splslevel.res.4$cor.value[1,1], 0.9637933)
})
