data(multidrug)
set.seed(42)
tune.spca.res <- tune.spca(X = multidrug$ABC.trans, ncomp = 2, 
                           nrepeat = 5, folds = 3, test.keepX = seq(5,15,5))
tune.spca.res
plot(tune.spca.res)
\dontrun{
## parallel processing using BiocParallel on repeats with 2 workers (cpus)
tune.spca.res <- tune.spca(X = multidrug$ABC.trans, ncomp = 2, 
                           nrepeat = 10, folds = 3, test.keepX = seq(5,15,5), 
                           ## use SnowParam() on Windows OS
                           BPARAM = MulticoreParam(workers = 2))
}
