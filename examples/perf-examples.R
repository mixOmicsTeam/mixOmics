## validation for objects of class 'pls' (regression)
# ----------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

# try tune the number of component to choose
# ---------------------
# first learn the full model
liver.pls <- pls(X, Y, ncomp = 5)

# with 5-fold cross validation: we use the same parameters as in model above
# but we perform cross validation to compute the MSEP, Q2 and R2 criteria
# ---------------------------
liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)

# see available criteria
names(liver.val$measures)
# see values for all repeats
liver.val$measures$Q2.total$values
# see summary over repeats
liver.val$measures$Q2.total$summary
# Q2 total should decrease until it reaches a threshold
liver.val$measures$Q2.total

# ncomp = 2 is enough
plot(liver.val, criterion  = 'Q2.total')

\dontrun{

# have a look at the other criteria
# ----------------------
# R2
plot(liver.val, criterion  = 'R2')
## correlation of components (see docs)
plot(liver.val, criterion  = 'cor.tpred')

# MSEP
plot(liver.val, criterion  = 'MSEP')
## validation for objects of class 'spls' (regression)
# ----------------------------------------
ncomp = 7
# first, learn the model on the whole data set
model.spls = spls(X, Y, ncomp = ncomp, mode = 'regression',
                  keepX = c(rep(10, ncomp)), keepY = c(rep(4,ncomp)))


# with leave-one-out cross validation
set.seed(45)
model.spls.val <- perf(model.spls, validation = "Mfold", folds = 5 )

#Q2 total
model.spls.val$measures$Q2$summary

# R2: we can see how the performance degrades when ncomp increases
plot(model.spls.val, criterion="R2")

## validation for objects of class 'splsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class

ncomp = 2

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
dist = "all", auc = TRUE)
error
error$auc

plot(error)

# parallel code
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
dist = "all", auc = TRUE, cpus =2)

# with 5 components and nrepeat=5, to get a $choice.ncomp
ncomp = 5
srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))

set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
dist = "all", nrepeat =5)
error$choice.ncomp

plot(error)


## validation for objects of class 'mint.splsda' (classification)
# ----------------------------------------

data(stemcells)
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, 
                  ncomp = 3, keepX = c(10, 5, 15),
                  study = stemcells$study)

out = perf(res, auc = TRUE)
out
plot(out)
out$auc
out$auc.study

## validation for objects of class 'sgccda' (classification)
# ----------------------------------------

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)

nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = 'full',
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,
scheme = "horst")

perf = perf(nutrimouse.sgccda)
perf
plot(perf)


# with 5 components and nrepeat=5 to get $choice.ncomp
nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = 'full',
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 5,
scheme = "horst")

perf = perf(nutrimouse.sgccda, folds = 5, nrepeat = 5)
perf
plot(perf)
perf$choice.ncomp
}