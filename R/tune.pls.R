tune.pls <-  function(X,
                      Y,
                      ncomp,
                      validation = c('Mfold', 'loo'),
                      nrepeat = 1,
                      folds,
                      mode = c('regression', 'canonical', 'classic'),
                      measure.tune = c("Q2.total", "RSS", "cor"),
                      BPPARAM = SerialParam(),
                      progressBar = FALSE,
                      ...
)
    {
    # TODO # tune.pls should be like before -- no cor and RSS (MAE, MSE, Bias R2)
    ## eval all args
    mget(names(formals()), sys.frame(sys.nframe()))

    ## run pls
    pls.res <- pls(X = X, Y = Y, ncomp = ncomp, mode = mode, all.outputs = FALSE, ...)
    ## evaluate performance
    result <- perf.mixo_pls(object = pls.res, 
                              validation = validation, 
                              folds = folds, 
                              progressBar = progressBar, 
                              nrepeat = nrepeat)
    meaure <- result[[measure.tune]]
    result$call <- match.call()
        
        return(result)
        
    }