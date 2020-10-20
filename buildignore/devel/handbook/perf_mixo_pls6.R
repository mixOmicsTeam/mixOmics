# changes to bypass the loop for the Q2

# Al: include check entries
# ! remove invariant as we cannot calculate RSS

perf.mixo_pls6 <- function(object,
                           validation = c("Mfold", "loo"),
                           folds = 10,
                           progressBar = FALSE,
                           ...)
  
{
  #------------------#
  #-- check entries --#
  
  # #-- check spls mode
  # if (object$mode == 'canonical')
  #   stop("'perf' is only available for (s)pls with modes: 'regression', 'invariant' or 'classic'.  Object has mode 'canonical'", call. = FALSE)
  # 
  # #-- validation
  # choices = c("Mfold", "loo")
  # validation = choices[pmatch(validation, choices)]
  # 
  # if (any(is.na(validation)) || length(validation) > 1)
  #   stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
  # 
  # #-- progressBar
  # if (!is.logical(progressBar))
  #   stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
  # 
  # #-- end checking --#
  # #------------------#
  
  
  #-- cross-validation approach  ---------------------------------------------#
  #---------------------------------------------------------------------------#
  
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
  
  if (any(is.na(X)) || any(is.na(Y)))
    stop("missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.", call. = FALSE)
  
  
  #-- tells which variables are selected in X and in Y --#
  if (is(object, "mixo_spls"))
  {
    keepX = object$keepX
    keepY = object$keepY
  } else {
    keepX = rep(ncol(X), ncomp)
    keepY = rep(ncol(Y), ncomp)
  }
  
  #-- define the folds --#
  if (validation == "Mfold")
  {
    if (is.list(folds))
    {
      
      if (length(folds) < 2 || length(folds) > n)
        stop("Invalid number of folds.", call. = FALSE)
      
      if (length(unlist(folds)) != n)
        stop("Invalid folds. The total number of samples in folds must be equal to ",
             n, ".", call. = FALSE)
      
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds. Repeated samples in folds.", call. = FALSE)
      
      M = length(folds)
    } else {
      if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n)
      {
        stop("Invalid number of folds.", call. = FALSE)
      } else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n))
      }
    }
  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
  }
  
  #-- set up a progress bar --#
  if (progressBar == TRUE)
  {
    pb = txtProgressBar(style = 3)
    nBar = 1
  }
  
  #-- initialize new objects --#
  if (mode == 'canonical'){
    RSS = rbind(rep(n - 1, p), matrix(nrow = ncomp, ncol = p))
    # RSS.indiv is the reconstructed matrix X
    #RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = p)})
    #RSS.indiv[[1]] = X
    press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = p)})
    PRESS.inside = Q2 = matrix(nrow = ncomp, ncol = p)
  }else{
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    # RSS.indiv is the reconstructed matrix Y
    #RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
    #RSS.indiv[[1]] = Y # KA changed
    press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = q)})
    PRESS.inside = Q2 = matrix(nrow = ncomp, ncol = q)
  }
  
  MSEP.mat = Ypred = array(0, c(n, q, ncomp))
  MSEP = R2 = matrix(nrow = ncomp, ncol = q)
  
  # to store the predicted components
  t.pred.cv = matrix(nrow = nrow(X), ncol = ncomp)
  u.pred.cv = matrix(nrow = nrow(X), ncol = ncomp)
  
  # to record feature stability 
  featuresX  = featuresY =  list()
  for(k in 1:ncomp){
    featuresX[[k]] = featuresY[[k]] = NA
  }
  
  # ====  loop on h = ncomp is only for the calculation of Q2 on each component
  for (h in 1:ncomp)
  {
    #-- initialising arguments --#
    tt = object$variates$X[, h]
    u = object$variates$Y[, h]
    b = object$loadings$Y[, h]
    #nx = p - keepX[h]
    #ny = q - keepY[h]
    
    # only used for matrices deflation across dimensions
    c = crossprod(X, tt)/drop(crossprod(tt))  #object$mat.c[, h]
    d = crossprod(Y, tt)/drop(crossprod(tt))  #object$mat.d[, h]
    e = crossprod(Y, u)/drop(crossprod(u))    
    
    # deflate matrices
    X = X - tt %*% t(c)
    
    #-- mode classic
    if (mode == "classic")
      Y = Y - tt %*% t(b)
    #-- mode regression
    if (mode == "regression")
      Y = Y - tt %*% t(d)
    #-- mode canonical 
    if (mode == "canonical")
      Y = Y - u %*% t(e)
    #-- mode invariant: Y is unchanged
    
    # update RSS for X/Y deflated
    if(mode == 'canonical'){  # based on X
      RSS[h + 1, ] =  colSums((X)^2)   # ==  colSums((X - tt %*% t(c))^2) if we had not deflated
    }else{ # regression, invariant, classic
      RSS[h + 1, ] = colSums((Y)^2)  # 
    }
    
  } # end h to calculate RSS   
  
  
  
  # ======== loop on i for cross validation ===================#
  for (i in 1:M)
  {
    if (progressBar == TRUE)
    {
      setTxtProgressBar(pb, nBar/(ncomp * M))
      nBar = nBar + 1
    }
    
    # initialise the train / test datasets
    omit = folds[[i]]
    X.train = object$X[-omit, , drop = FALSE]
    Y.train = object$Y[-omit, , drop = FALSE]
    X.test = object$X[omit, , drop = FALSE]
    Y.test = object$Y[omit, , drop = FALSE]
    
    # New loop to calculate prediction
    for (h in 1:ncomp)
    { 
      #-- for MSEP and R2 criteria, no loop on the component as we do a spls with ncomp
      ##if (h == 1)
      #{
      #nzv = (apply(X.train, 2, var) > .Machine$double.eps) # removed in v6.0.0 so that MSEP, R2 and Q2 are obtained with the same data
      # re-added in >6.1.3 to remove constant variables
      nzv = (apply(X.train, 2, var) > .Machine$double.eps)
      
      # creating a keepX.temp that can change for each fold, depending on nzv
      keepX.temp = keepX
      if(any(keepX.temp > sum(nzv)))
        keepX.temp[which(keepX.temp>sum(nzv))] = sum(nzv)
      
      # here h = 1 because we deflate at each step then extract the vectors for each h comp
      spls.res = mixOmics::spls(X.train[,nzv], Y.train, ncomp = 1, mode = mode, max.iter = max.iter, tol = tol, 
                                keepX = keepX.temp, keepY = keepY, near.zero.var = FALSE, scale = scale)
      Y.hat = mixOmics:::predict.mixo_spls(spls.res, X.test[,nzv, drop = FALSE])$predict
      
      # added the stop msg
      if(sum(is.na(Y.hat))>0) stop('Predicted Y values include NA')  
      
      # replaced h by 1; Y.hat is the prediction of the test samples for all q variable in comp h = 1
      Ypred[omit, , h] = Y.hat[, , 1]
      MSEP.mat[omit, , h] = (Y.test - Y.hat[, , 1])^2
      
      
      # Q2 criterion: buidling directly from spls object
      u.cv = spls.res$variates$Y[, 1]
      t.cv = spls.res$variates$X[, 1]
      a.cv = spls.res$loadings$X[, 1]
      b.cv = spls.res$loadings$Y[, 1]
      
      # reg coefficients:
      c.cv = crossprod(X.train, u.cv) / drop(crossprod(u.cv)) 
      d.cv = crossprod(Y.train, t.cv) / drop(crossprod(t.cv)) # d.cv \neq to b.cv as d.cv is normed wrt to t.cv
      e.cv = crossprod(Y.train, u.cv) / drop(crossprod(u.cv)) 
      
      # calculate predicted components and store
      t.pred = c(X.test %*% a.cv)
      t.pred.cv[omit,h] = t.pred    # needed for tuning
      b.pred = crossprod(Y.test, t.pred)
      b.pred.cv = b.pred/ drop(sqrt(crossprod(b.pred)))
      u.pred.cv[omit,h] = Y.test %*% b.cv  # needed for tuning, changed instead of b.pred.cv
      
      # predicted reg coeff, could be removed
      e.pred.cv = crossprod(as.matrix(Y.test), Y.test %*% b.pred.cv) / drop(crossprod(Y.test %*% b.pred))
      d.pred.cv = crossprod(as.matrix(Y.test), t.pred) / drop(crossprod(t.pred)) 
      
      # deflate matrices X
      X.train = X.train - t.cv %*% t(c.cv)
      X.test = X.test - t.pred %*% t(c.cv)
      # deflate matrices X      
      #-- mode classic
      if (mode == "classic"){
        Y.train = Y.train - t.cv %*% t(b.cv)  # could be pred on b
        Y.test = Y.test - t.pred %*% t(b.cv)
      }
      #-- mode regression
      if (mode == "regression"){
        Y.train = Y.train - t.cv %*% t(d.cv) # could be pred d.pred.cv? does not decrease enough
        Y.test = Y.test - Y.hat[, , 1]   # == Y.test - t.pred %*% t(d.cv) 
      }
      #-- mode canonical  ## KA added
      if (mode == "canonical"){
        Y.train = Y.train - u.cv %*% t(e.cv)  # could be pred on e
        Y.test = Y.test - (Y.test %*% b.cv) %*% t(e.cv)  # here u.pred = Y.test %*% b.cv (b.pred.cv gives similar results)
      }
      #-- mode invariant: Y is unchanged
      
      # calculate predicted matrix X.hat or Y.hat based on X.test
      if(mode == 'canonical'){  # Xa c' = t c'
        #X.hat.cv = t.pred %*% t(c.cv), calculated earlier
        press.mat[[h]][omit, ] = X.test        # == X.test - X.hat.cv
      }else{ #  if(mode == 'regression'){  # Xa d' = t d'
        #Y.hat.cv = t.pred %*% t(d.cv), calculated earlier
        press.mat[[h]][omit, ] = Y.test        # == Y.test - Y.hat.cv
      }  
      
      
      # Record selected features in each set
      if (is(object,"mixo_spls"))
      {
        featuresX[[h]] = c(unlist(featuresX[[h]]), selectVar(spls.res, comp = 1)$X$name)
        featuresY[[h]] = c(unlist(featuresY[[h]]), selectVar(spls.res, comp = 1)$Y$name)
      }
      
    } #  end loop on h ncomp
  } # end i (cross validation)
  
  
  
  # store results for each comp
  for (h in 1:ncomp){
    #-- compute the Q2 criterion --#
    # norm is equivalent to summing here the squared press values:
    PRESS.inside[h, ] = apply(press.mat[[h]], 2, function(x){norm(x, type = "2")^2})
    
    if(mode != 'canonical'){
    Q2[h, ] = 1 - PRESS.inside[h, ] / RSS[h, ]
    MSEP[h, ] = apply(as.matrix(MSEP.mat[, , h]), 2, mean)
    R2[h, ] = (diag(cor(object$Y, Ypred[, , h])))^2
    } # if mode == canonical, do not output
  }
  
  
  if (progressBar == TRUE) cat('\n')
  
  
  #-- output -----------------------------------------------------------------#
  #---------------------------------------------------------------------------#
  Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]),
                    nrow = 1, ncol = ncomp,
                    dimnames = list("Q2.total", paste0(1:ncomp, " comp")))
  
  # set up dimnames and outputs
  result = list()
  
  if(mode != 'canonical'){
  rownames(MSEP) = rownames(R2) = rownames(Q2) = paste0(1:ncomp, " comp")
  colnames(MSEP) = colnames(R2) = colnames(Q2) = object$names$colnames$Y
  
  result$MSEP = t(MSEP)
  #result$MSEP.mat = MSEP.mat  
  result$R2 = t(R2)
  result$Q2 = t(Q2)  # remove this output when canonical mode?
  }
  
  

  result$Q2.total =  Q2.total
  result$RSS = RSS
  result$PRESS = PRESS.inside
  # result$press.mat = press.mat
  #result$RSS.indiv = RSS.indiv
  result$d.cv = d.cv  # KA added  
  result$b.cv = b.cv  # KA added 
  result$c.cv = c.cv  # KA added 
  result$u.cv = u.cv  # KA added 
  result$a.cv = a.cv  # KA added 
  result$t.pred.cv = t.pred.cv  # needed for tuning
  result$u.pred.cv = u.pred.cv  # needed for tuning

  
  #---- extract stability of features -----#
  if (is(object, "mixo_spls"))
  {
    list.features.X = list()
    list.features.Y = list()
    
    for(k in 1:ncomp)
    {
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.features.X[[k]] = sort(table(as.factor(featuresX[[k]][-remove.naX])) / M, decreasing = TRUE)
      list.features.Y[[k]] = sort(table(as.factor(featuresY[[k]][-remove.naY])) / M, decreasing = TRUE)
      
    }
    names(list.features.X)  = names(list.features.Y) = paste0('comp', 1:ncomp)
    
    # features
    result$features$stable.X = list.features.X
    result$features$stable.Y = list.features.Y
  }
  
  #--- class
  if (is(object,"mixo_spls"))
  {
    method = "spls.mthd"
  } else if (is(object, "mixo_pls")) {
    method = "pls.mthd"
  } else {
    warning("Something happened. Please contact us.")
  }
  class(result) = c("perf",paste(c("perf", method), collapse ="."))
  result$call = match.call()
  
  return(invisible(result))
}
