
source('perf_mixo_pls6.R')

# ! type.tune shoud be null for a PLS model

tune_spls2_repeat = function(X, Y, list.keepX, list.keepY, ncomp, nrepeat, mode, type.tune, pls.model){
  out = list()
  
  if(isFALSE(pls.model)){
    cor.tpred = cor.upred = RSS.tpred = RSS.upred = array(dim = c(length(list.keepX), length(list.keepY), nrepeat),
                                                          dimnames = list(paste0('keepX', list.keepX), 
                                                                          paste0('keepY', list.keepY),
                                                                          paste0('repeat', 1:nrepeat)))
    best.keepX = best.keepY = NULL
  }else{
    type.tune = NULL
    cor.tpred = cor.upred = RSS.tpred = RSS.upred = matrix(nrow = ncomp, ncol = nrepeat, 
                                                           dimnames = list(paste0('comp', 1:ncomp), paste0('repeat', 1:nrepeat)))
    Q2.tot.ave = matrix(nrow = ncomp, ncol = nrepeat,
                        dimnames = list(paste0('comp', 1:ncomp), paste0('repeat', 1:nrepeat)))
    
  }
  cor.pred = RSS.pred = list()
  cov.pred = list()
  
  
  
  ## for a PLS only to extract Q2.total (or anything else)
  if(isTRUE(pls.model)){
      for(k in 1:nrepeat){
      cat('repeat', k, '\n')
      
      # run a full PLS model up to the last comp
      pls.res = spls(X = X, Y = Y, 
                     keepX = c(rep(ncol(X), ncomp)), 
                     keepY = c(rep(ncol(Y), ncomp)), 
                     ncomp = ncomp, mode = mode)
      # fold CV
      res.perf = perf.mixo_pls6(pls.res, validation = 'Mfold', folds = 10)
      
      # extract Q2.total for a PLS, we could extract other outputs such as R2, MSEP etc (only valid for regression)
      Q2.tot.ave[, k] = res.perf$Q2.total 
      
      # extract the predicted components per dimension, take abs value
      cor.tpred[, k] = diag(abs(cor(res.perf$t.pred.cv, pls.res$variates$X)))
      cor.upred[, k] = diag(abs(cor(res.perf$u.pred.cv, pls.res$variates$Y)))
      
      # RSS: no abs values here
      RSS.tpred[, k] = apply((res.perf$t.pred.cv - pls.res$variates$X)^2, 2, sum)/(nrow(X) -1)
      RSS.upred[, k] = apply((res.perf$u.pred.cv - pls.res$variates$Y)^2, 2, sum)/(nrow(X) -1)
    } #end repeat       
    
    # # calculate mean across repeats
    cor.pred$u = apply(cor.upred, 1, mean) 
    cor.pred$t = apply(cor.tpred, 1, mean)
    RSS.pred$u = apply(RSS.upred, 1, mean)
    RSS.pred$t = apply(RSS.tpred, 1, mean)
    
  }else{ # if sPLS model 
    for (comp in 1:ncomp){
      cat('Comp', comp, '\n')
      for(k in 1:nrepeat){
        cat('repeat', k, '\n')
        for(keepX in 1:length(list.keepX)){
          #cat('KeepX', list.keepX[keepX], '\n')
          for(keepY in 1:length(list.keepY)){
            #cat('KeepY', list.keepY[keepY], '\n')
            
            # sPLS model, updated with the best keepX
            spls.res = spls(X = X, Y = Y, 
                            keepX = c(best.keepX, list.keepX[keepX]), 
                            keepY = c(best.keepY, list.keepY[keepY]), 
                            ncomp = comp, mode = mode)
            # fold CV
            res.perf = perf.mixo_pls6(spls.res, validation = 'Mfold', folds = 10)
            
            # extract the predicted components: 
            if(type.tune == 'cor' ){
              cor.tpred[keepX, keepY, k] = cor(res.perf$t.pred.cv[,comp], spls.res$variates$X[, comp])
              cor.upred[keepX, keepY,k] = cor(res.perf$u.pred.cv[,comp], spls.res$variates$Y[, comp])
            }
            if(type.tune == 'RSS'){
              # RSS: no abs values here
              RSS.tpred[keepX, keepY, k] = sum((res.perf$t.pred.cv[,comp] - spls.res$variates$X[, comp])^2)/(nrow(X) -1) 
              RSS.upred[keepX, keepY, k] = sum((res.perf$u.pred.cv[,comp] - spls.res$variates$Y[, comp])^2)/(nrow(X) -1)
            }
            # covariance between predicted variates
            ##cov.variate.pred[keepX, keepY, k] = cov(res.perf$t.pred.cv[,comp], res.perf$u.pred.cv[,comp])
          } # end keepY
        } #end keepX
      } # end repeat
      cat('\t')
      
      # # calculate mean across repeats
      cor.pred$u[[comp]] = apply(cor.upred, c(1,2), mean)  #mean(cor.upred[keepX, keepY,])
      cor.pred$t[[comp]] = apply(cor.tpred, c(1,2), mean)  #mean(cor.tpred[keepX, keepY,])
      RSS.pred$u[[comp]] = apply(RSS.upred, c(1,2), mean)  #mean(RSS.upred[keepX, keepY,])
      RSS.pred$t[[comp]] = apply(RSS.tpred, c(1,2), mean)  #mean(RSS.tpred[keepX, keepY,])
      
      # choose the best keepX and keepY based on type.tune
      if(mode != 'canonical'){  #regression, invariant, classic
        # define best keepX and keepY based on u
        if(type.tune == 'cor'){
          cor.component = cor.pred$u[[comp]]
          index = which(cor.component == max(cor.component), arr.ind = TRUE)
        }else{ # if type.tune = 'RSS'
          RSS.component = RSS.pred$u[[comp]]
          index = which(RSS.component == min(RSS.component), arr.ind = TRUE)
        }
        best.keepX = c(best.keepX, list.keepX[index[1,1]])
        best.keepY = c(best.keepY, list.keepY[index[1,2]])
        
      }else{  # mode = 'canonical'
        if(type.tune == 'cor'){
          cor.component.t = cor.pred$t[[comp]]
          cor.component.u = cor.pred$u[[comp]]
          index.t = which(cor.component.t == max(cor.component.t), arr.ind = TRUE)
          index.u = which(cor.component.u == max(cor.component.u), arr.ind = TRUE)
        }else{ # if type.tune = 'RSS'
          RSS.component.t = RSS.pred$t[[comp]]
          RSS.component.u = RSS.pred$u[[comp]]
          index.t = which(RSS.component.t == min(RSS.component.t), arr.ind = TRUE)
          index.u = which(RSS.component.u == min(RSS.component.u), arr.ind = TRUE)
        }
        best.keepX = c(best.keepX, list.keepX[index.t[1,1]])
        best.keepY = c(best.keepY, list.keepY[index.u[1,2]])
      } # canonical
    } # end comp
    
    out$best.keepX = best.keepX
    out$best.keepY = best.keepY
  } # end sPLS
  
  
  
  
  # 
  # # # summary of results
  # summary = list()
  # for(k in 1:ncomp){
  #   RSS.upred = RSS.pred$u[[k]]
  #   RSS.tpred = RSS.pred$t[[k]]
  # 
  #   cor.upred = cor.pred$u[[k]]
  #   cor.tpred = cor.pred$t[[k]]
  # 
  #   # u.pred
  #   index = which(cor.upred == max(cor.upred), arr.ind = TRUE)
  #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
  #   #cor.upred[index]
  #   summary$u[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], cor.upred[index])
  # 
  #   # t.pred
  #   index = which(cor.tpred == max(cor.tpred), arr.ind = TRUE)
  #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
  #   #cor.tpred[index]
  #   summary$t[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], cor.tpred[index])
  # 
  #   # RSS.u: same as upred
  #   index = which(RSS.upred == min(RSS.upred), arr.ind = TRUE)
  #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
  #   #RSS.upred[index]
  #   summary$RSSu[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], RSS.upred[index])
  # 
  # 
  #   # RSS.t: different from tpred?
  #   index = which(RSS.tpred == min(RSS.tpred), arr.ind = TRUE)
  #   #list.keepX[index[1,1]]; list.keepY[index[1, 2]];
  #   #RSS.tpred[index]
  #   summary$RSSt[[k]] = c(list.keepX[index[1,1]], list.keepY[index[1, 2]], RSS.tpred[index])
  # } # end ncomp
  #   out$summary = summary
  
  if(isTRUE(pls.model)){
    out$cor.pred = cor.pred
    out$RSS.pred = RSS.pred
    out$Q2.tot.ave = apply(Q2.tot.ave, 1, mean)
  }else{
    if(type.tune == 'cor') out$cor.pred = cor.pred
    if(type.tune == 'RSS') out$RSS.pred = RSS.pred
  }
  
  return(out)
} 