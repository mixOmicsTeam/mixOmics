data("nutrimouse")
nrepeat <- 5
tune.spca.res <- tune.spca(
    X = nutrimouse$lipid,
    ncomp = 2,
    nrepeat = nrepeat,
    folds = 3,
    test.keepX = seq(5, 15, 5),
    seed = 42
)
tune.spca.res
plot(tune.spca.res)
\dontrun{
## parallel processing using BiocParallel on repeats with more workers (cpus)
# Check if the environment variable exists (during R CMD check) and limit cores accordingly
max_cores <- if (Sys.getenv("_R_CHECK_LIMIT_CORES_") != "") 2 else parallel::detectCores() - 1
# Setup the parallel backend with the appropriate number of workers
BPPARAM <- BiocParallel::MulticoreParam(workers = max_cores)
tune.spca.res <- tune.spca(
    X = nutrimouse$lipid,
    ncomp = 2,
    nrepeat = nrepeat,
    folds = 3,
    test.keepX = seq(5, 15, 5),
    BPPARAM = BPPARAM
)
plot(tune.spca.res)
}
