#' Tune number of selected variables for spca
#'
#' @inheritParams spca
#' @inheritParams tune.splsda
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating the type
#'   of parallelisation. See examples.
#' @importFrom BiocParallel SerialParam bplapply
#' @return A \code{tune.spca} object containing: \describe{ 
#' \item{call}{ The
#'   function call}
#' \item{choice.keepX}{The selected number of components on
#'   each component} 
#' \item{cor.comp}{The correlations between the components
#'   from the cross-validated studies and those from the study which used all of
#'   the data in training.} }
#' @export
#'
#' @example ./examples/tune.spca-examples.R
tune.spca <- function(X, ncomp=2, nrepeat=3, folds, test.keepX, center = TRUE, scale = TRUE, 
                      BPPARAM = SerialParam())
{
    if (nrepeat < 3)
    {
        stop("'nrepeat' must be >= 3")
    }
    ## optimal keepX for all components
    keepX.opt <- NULL
    ## a list of cor.df for each component
    cor.df.list <- .name_list(char = paste0('comp', seq_len(ncomp)))
    # cor.df.list <- lapply(cor.df.list, function(x) cor.df)
    
    if (nrow(X) < 30) {
        cat("\n Low sample size. Using Leave-One-Out CV without repeat")
        nrepeat <- 1
        folds <- nrow(X)
    }
    all.keepX <- test.keepX
    names(all.keepX) <- paste0('keepX_', all.keepX)
    ## ------ component loop
    for(ncomp in seq_len(ncomp)) {
        iter_keepX <- function(keepX.value) {
            ## ------ repeated cv
            repeat_cv_j <- function(j) {
                repeat.j.folds = split(sample(seq_len(nrow(X))),seq_len(folds))
                ## ------ mean cor for CV
                cor.pred = sapply(repeat.j.folds, function(test.fold.inds){
                    ## split data to train/test
                    X.train = X[-test.fold.inds,,drop=FALSE]
                    X.test = X[test.fold.inds,,drop=FALSE]
                    X.test = scale(X.test, center = center, scale = scale)
                    
                    # ---- run sPCA 
                    suppressWarnings({
                        ## train
                        spca.train = mixOmics::spca(X.train, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = scale)
                        # full data
                        spca.full = mixOmics::spca(X, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = scale)
                    })
                    # ---- deflation on X with only the fold left out
                    for(k in seq_len(ncomp)){ 
                        # loop to calculate deflated matrix and predicted comp
                        # calculate the predicted comp on the fold left out
                        # calculate reg coeff, then deflate
                        if(k == 1){
                            t.comp.sub = X.test %*% spca.train$loadings$X[,k]
                        } else{
                            # calculate deflation beyond comp 1
                            # recalculate the loading vector (here c.sub) on the test set (perhaps we could do this instead on the training set by extracting from spca.train$loadings$X[,k]?)
                            c.sub = crossprod(X.test, t.comp.sub) / drop(crossprod(t.comp.sub)) 
                            X.test = X.test - t.comp.sub %*% t(c.sub) 
                            # update predicted comp based on deflated matrix
                            t.comp.sub = X.test %*% spca.train$loadings$X[,k]
                        }
                    }
                    
                    # calculate predicted component and compare with component from sPCA on full data on the left out set
                    # cor with the component on the full data, abs value
                    cor.comp = abs(cor(t.comp.sub, spca.full$variates$X[,ncomp][test.fold.inds], use = 'pairwise.complete.obs'))
                    # need to flip the sign in one of the comp to calculate the RSS
                    #RSS.comp = sum((c(sign(cor.comp))*(t.comp.sub) - spca.full$variates$X[,comp][i])^2)
                    # also look at correlation between loading vectors
                    # cor.loading = abs(cor(spca.full$loadings$X[,comp], spca.train$loadings$X[,comp]))
                    # get the cor(pred component, comp) per fold for a given comp
                    return(cor.comp)
                })  # end CV, get the cor(pred component, comp) per fold for a given comp
                
                # output is correlation between reconstructed components cor.comp on the left out samples and if we had those samples in PCA,
                # average correlations across folds for the repeat
                return(mean(cor.pred, na.rm = TRUE))
            }
            out <-  lapply(seq_len(nrepeat), FUN=repeat_cv_j)
            return(unlist(out))
        }
        test.keepX.cors <- bplapply(X=all.keepX, FUN=iter_keepX, BPPARAM = BPPARAM)
        test.keepX.cors <- data.frame(test.keepX.cors)
        test.keepX.cors <- t(test.keepX.cors)
        colnames(test.keepX.cors) <- paste0('repeat_', seq_len(nrepeat))
        cor.df.list[[ncomp]] <- test.keepX.cors
        # browser()
        ## use a one-sided t.test using repeat correlations to assess if addition of keepX improved the correlation
        ## and get the index of optimum keepX
        t.test.df <- data.frame(t(cor.df.list[[ncomp]]))
        keepX.opt.comp.ind <-  t.test.process(t.test.df, alpha = 0.1, alternative = 'less')
        keepX.opt.comp <- all.keepX[keepX.opt.comp.ind]
        ## update keepX.optimum for next comp
        keepX.opt <- c(keepX.opt, keepX.opt.comp)
    }
    choice.keepX <- keepX.opt
    names(choice.keepX) <- paste0('comp', seq_len(ncomp))
    
    ## get the mean and sd values
    cor.comp = mapply(df=cor.df.list, opt = choice.keepX, FUN = function(df, opt){
        out <- data.frame(keepX = all.keepX, cor.mean = apply(df, 1, mean), cor.sd = apply(df, 1, sd))
        out$opt.keepX <- NA
        out[which(out$keepX == opt),]$opt.keepX <- TRUE
        return(out)
    }, SIMPLIFY = FALSE)
    # evaluate all for output except X to save memory
    mc <- mget(names(formals())[-1], sys.frame(sys.nframe()))
    mc <- as.call(c(as.list(match.call())[1:2], mc))
    
    result <- list(
        call = mc,
        choice.keepX = choice.keepX,
        cor.comp = cor.comp)
    class(result) <- 'tune.spca'
    return(result)
}
