
\dontrun{
## validation for objects of class 'splsda'

data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)
out = tune.splsda(X, Y, ncomp = 3, nrepeat = 5, logratio = "none",
test.keepX = c(5, 10, 15), folds = 10, dist = "max.dist",
progressBar = TRUE)


plot(out, sd=TRUE)
}
\dontrun{
## validation for objects of class 'mint.splsda'

data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study

out = tune(method="mint.splsda", X=data,Y=type.id, ncomp=2, study=exp, test.keepX=seq(1,10,1))
out$choice.keepX

plot(out)

## validation for objects of class 'mint.splsda'

data("breast.TCGA")
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
data = with(breast.TCGA$data.train, list(mrna = mrna, 
            mirna = mirna,
            protein = protein,
            Y = subtype))
# set number of component per data set
ncomp = 5

# Tuning the first two components
# -------------

# definition of the keepX value to be tested for each block mRNA miRNA and protein
# names of test.keepX must match the names of 'data'
test.keepX = list(mrna = seq(10,40,20), mirna = seq(10,30,10), protein = seq(1,10,5))

# the following may take some time to run, note that for through tuning
# nrepeat should be > 1
tune = tune.block.splsda(X = data, indY = 4,
ncomp = ncomp, test.keepX = test.keepX, design = 'full', nrepeat = 3)

tune$choice.ncomp
tune$choice.keepX

plot(tune)
## --- spls model
data(nutrimouse)
X <- nutrimouse$gene  
Y <- nutrimouse$lipid
list.keepX <- c(2:10, 15, 20)
# tuning based on correlations
set.seed(30)
## tune X only
tune.spls.cor.X <- tune.spls(X, Y, ncomp = 3,
                           test.keepX = list.keepX,
                           validation = "Mfold", folds = 5,
                           nrepeat = 3, progressBar = FALSE,
                           measure = 'cor')

plot(tune.spls.cor.X)
plot(tune.spls.cor.X, measure = 'RSS')

## tune Y only
tune.spls.cor.Y <- tune.spls(X, Y, ncomp = 3,
                           test.keepY = list.keepX,
                           validation = "Mfold", folds = 5,
                           nrepeat = 3, progressBar = FALSE,
                           measure = 'cor')

plot(tune.spls.cor.Y)
plot(tune.spls.cor.Y, sd = FALSE)
plot(tune.spls.cor.Y, measure = 'RSS')

## tune Y and X
tune.spls.cor.XY <- tune.spls(X, Y, ncomp = 3,
                             test.keepY = c(8, 15, 20),
                             test.keepX = c(8, 15, 20),
                             validation = "Mfold", folds = 5,
                             nrepeat = 3, progressBar = FALSE,
                             measure = 'cor')

plot(tune.spls.cor.XY)
## show RSS
plot(tune.spls.cor.XY, measure = 'RSS')
## customise point sizes
plot(tune.spls.cor.XY, size.range = c(6,12))
}
