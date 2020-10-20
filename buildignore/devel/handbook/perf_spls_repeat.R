
source('perf_mixo_pls6.R')

# ! be mindful on how to implement this for multilevel?

perf_spls_repeat = function(object, validation, folds, nrepeat, measure){
  
  #-- initialising arguments --#
  # these are the centered and scaled matrices output from pls, we remove $nzv if needed
  if (length(object$nzv$Position)>0)
  {
    X = object$X[, -object$nzv$Position]
  } else {
    X = object$X
  }
  Y = object$Y
  
  scale = object$scale
  tol = object$tol
  max.iter = object$max.iter
  mode = object$mode
  ncomp = object$ncomp
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  
  keepX = object$keepX
  keepY = object$keepY
  
  
  out = list()
  
    measure = NULL
    cor.tpred = cor.upred = RSS.tpred = RSS.upred = matrix(nrow = ncomp, ncol = nrepeat, 
                                                           dimnames = list(paste0('comp', 1:ncomp), paste0('repeat', 1:nrepeat)))
    Q2.tot.ave = matrix(nrow = ncomp, ncol = nrepeat,
                        dimnames = list(paste0('comp', 1:ncomp), paste0('repeat', 1:nrepeat)))

  cor.pred = RSS.pred = list()
 
    for(k in 1:nrepeat){
      cat('repeat', k, '\n')
      # fold CV
      res.perf = perf.mixo_pls6(object, validation = validation, folds = folds)
      
      # extract Q2.total, we could extract other outputs such as R2, MSEP etc (only valid for regression)
      Q2.tot.ave[, k] = res.perf$Q2.total 
      
      # extract the predicted components per dimension, take abs value
      cor.tpred[, k] = diag(abs(cor(res.perf$t.pred.cv, object$variates$X)))
      cor.upred[, k] = diag(abs(cor(res.perf$u.pred.cv, object$variates$Y)))
      
      # RSS: no abs values here
      RSS.tpred[, k] = apply((res.perf$t.pred.cv - object$variates$X)^2, 2, sum)/(nrow(X) -1)
      RSS.upred[, k] = apply((res.perf$u.pred.cv - object$variates$Y)^2, 2, sum)/(nrow(X) -1)
    } #end repeat       
    
    # # calculate mean across repeats
    cor.pred$u = apply(cor.upred, 1, mean) 
    cor.pred$t = apply(cor.tpred, 1, mean)
    RSS.pred$u = apply(RSS.upred, 1, mean)
    RSS.pred$t = apply(RSS.tpred, 1, mean)
    
  # }else{ # if sPLS model 
  #   for(k in 1:nrepeat){
  #     cat('repeat', k, '\n')
  #     
  #     # fold CV
  #     res.perf = perf.mixo_pls6(object, validation = 'Mfold', folds = 10)
  #     
  #     # extract the predicted components: 
  #     if(measure == 'cor' ){
  #       cor.tpred[keepX, keepY, k] = cor(res.perf$t.pred.cv[,comp], object$variates$X[, comp])
  #       cor.upred[keepX, keepY,k] = cor(res.perf$u.pred.cv[,comp], object$variates$Y[, comp])
  #     }
  #     if(measure == 'RSS'){
  #       # RSS: no abs values here
  #       RSS.tpred[keepX, keepY, k] = sum((res.perf$t.pred.cv[,comp] - object$variates$X[, comp])^2)/(nrow(X) -1) 
  #       RSS.upred[keepX, keepY, k] = sum((res.perf$u.pred.cv[,comp] - object$variates$Y[, comp])^2)/(nrow(X) -1)
  #     }
  #   } # end repeat
  #   cat('\t')
  #   
  #   # # calculate mean across repeats
  #   cor.pred$u[[comp]] = apply(cor.upred, c(1,2), mean)  #mean(cor.upred[keepX, keepY,])
  #   cor.pred$t[[comp]] = apply(cor.tpred, c(1,2), mean)  #mean(cor.tpred[keepX, keepY,])
  #   RSS.pred$u[[comp]] = apply(RSS.upred, c(1,2), mean)  #mean(RSS.upred[keepX, keepY,])
  #   RSS.pred$t[[comp]] = apply(RSS.tpred, c(1,2), mean)  #mean(RSS.tpred[keepX, keepY,])
  #   
  # } # end sPLS
  
  
  # if(isTRUE(pls.model)){
    out$cor.pred = cor.pred
    out$RSS.pred = RSS.pred
    out$Q2.tot.ave = apply(Q2.tot.ave, 1, mean)
  # }else{
  #   if(measure == 'cor') out$cor.pred = cor.pred
  #   if(measure == 'RSS') out$RSS.pred = RSS.pred
  # }
  
  return(out)
} 