data(multidrug)
tune.spca_res <- tune.spca(X = multidrug$ABC.trans, ncomp = 2, nrepeat = 3, kfold = 3, test.keepX = seq(5,35,5))
tune.spca_res
