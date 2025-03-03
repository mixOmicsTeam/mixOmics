## PLS-DA example

data(liver.toxicity) # rats gex and clinical measurements/treatments
unique(liver.toxicity$treatment$Treatment.Group) # 16 groups
length(liver.toxicity$treatment$Treatment.Group) # 64 samples

plsda.res <- plsda(liver.toxicity$gene, liver.toxicity$treatment$Treatment.Group, ncomp = 2)

performance <- perf.assess(plsda.res, 
                        # to make sure each fold has all classes represented
                           validation = "Mfold", folds = 3, nrepeat = 10, 
                           seed = 12) # for reproducibility, remove for analysis

performance$error.rate$BER

## sPLS-DA example

splsda.res <- splsda(liver.toxicity$gene, liver.toxicity$treatment$Treatment.Group,
                     keepX = c(25, 25), ncomp = 2)

performance <- perf.assess(splsda.res,
                           validation = "Mfold", folds = 3, nrepeat = 10, 
                           seed = 12)

performance$error.rate$BER # can see slight improvement in error rate over PLS-DA example

## PLS example

ncol(liver.toxicity$clinic) # 10 Y variables as output of PLS model

pls.res <- pls(liver.toxicity$gene, liver.toxicity$clinic, ncomp = 2)

performance <- perf.assess(pls.res, 
                           validation = "Mfold", folds = 3, nrepeat = 10,
                           seed = 12)

# see Q2 which gives indication of predictive ability for each of the 10 Y outputs
performance$measures$Q2$summary 

## sPLS example

spls.res <- spls(liver.toxicity$gene, liver.toxicity$clinic, ncomp = 2, keepX = c(50, 50))

performance <- perf.assess(spls.res, 
                           validation = "Mfold", folds = 3, nrepeat = 10, 
                           seed = 12)

# see Q2 which gives indication of predictive ability for each of the 10 Y outputs
performance$measures$Q2$summary


## block PLS-DA example

data("breast.TCGA")
mrna <- breast.TCGA$data.train$mrna
mirna <- breast.TCGA$data.train$mirna
data <- list(mrna = mrna, mirna = mirna)
design <- matrix(1, 
   ncol = length(data), 
   nrow = length(data), 
   dimnames = list(names(data), names(data)))
diag(design) <-  0

block.plsda.res <- block.plsda(X = data, Y = breast.TCGA$data.train$subtype, 
                               ncomp = 2, design = design)

performance <- perf.assess(block.plsda.res)

performance$error.rate.per.class # error rate per class per distance metric

## block sPLS-DA example

block.splsda.res <- block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                 ncomp = 2, design = design,
                                 keepX = list(mrna = c(8,8), mirna = c(8,8)))

performance <- perf.assess(block.splsda.res)

performance$error.rate.per.class

## MINT PLS-DA example

data("stemcells")

mint.plsda.res <- mint.plsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3,
                             study = stemcells$study)

performance <- perf.assess(mint.plsda.res)

performance$global.error$BER # global error per distance metric

## MINT sPLS-DA example

mint.splsda.res <- mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, 
                               keepX = c(10, 5, 15), study = stemcells$study)

performance <- perf.assess(mint.splsda.res)

performance$global.error$BER # error slightly higher in this sparse model verses non-sparse