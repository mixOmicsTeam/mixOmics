# to Al: load the data, then go to line 119 (which has been commented currently) and onwards. This is the code
# used in the handbook, for section 10.5.4
# anything that called 'tune_spls2_repeat' and 'perf_spls_repeat' is what you should run.
## ----------------------------------------------------------------------------------------
devtools::load_all()
# data(liver.toxicity)
# X <- liver.toxicity$gene
# Y <- liver.toxicity$clinic
# 
# 
# ## ----------------------------------------------------------------------------------------
# head(data.frame(rownames(X), rownames(Y)))
# 
# 
# ## ---- eval = FALSE, fig.show = 'hide', message = FALSE-----------------------------------
# ## pls.result <- pls(X, Y)     # 1 Run the method
# ## plotIndiv(pls.result)       # 2 Plot the samples
# ## plotVar(pls.result)         # 3 Plot the variables
# 
# 
# ## ---- eval = FALSE, fig.show = 'hide', message = FALSE-----------------------------------
# ## spls.result <- spls(X, Y, keepX = c(10, 20), keepY = c(3, 2), ncomp = 2)
# ## plotIndiv(spls.result)
# ## plotVar(spls.result)
# 
# 
# ## ----------------------------------------------------------------------------------------
# y = liver.toxicity$clinic[, "ALB.g.dL."]
# 
# 
# ## ----pls1-Q2, fig.cap='(ref:pls1-Q2)'----------------------------------------------------
# pls1.liver.for.tune = pls(X = X, Y = y, ncomp = 10, mode = 'regression')
# set.seed(33)  # for reproducibility with this handbook, remove otherwise
# Q2.pls1.liver = perf(pls1.liver.for.tune, validation = 'Mfold', folds = 10)
# plot(Q2.pls1.liver,criterion = 'Q2')
# 
# 
# ## ----spls1-MAE, fig.cap='(ref:spls1-MAE)'------------------------------------------------
# # set up grid of values: 
# list.keepX = c(5:10, seq(15, 50, 5))     
# # list.keepX  # to have a look at the keepX grid
# set.seed(33)  # for reproducibility with this handbook, remove otherwise
# tune.spls1.MAE = tune.spls(X, y, ncomp= 2, 
#                           test.keepX = list.keepX, 
#                           validation = 'Mfold', folds = 10,
#                           nrepeat = 10, progressBar = FALSE, 
#                           measure = 'MAE')
# plot(tune.spls1.MAE)
# 
# 
# ## ----------------------------------------------------------------------------------------
# choice.ncomp = tune.spls1.MAE$choice.ncomp$ncomp
# # optimal number of variables to select in X based on the MAE criterion
# choice.keepX = tune.spls1.MAE$choice.keepX[1:choice.ncomp]  # we stop at choice.ncomp
# 
# choice.ncomp
# choice.keepX
# 
# 
# ## ----------------------------------------------------------------------------------------
# spls1.liver <- spls(X, y, ncomp = choice.ncomp, keepX = choice.keepX, mode = "regression")
# 
# 
# ## ---- eval = FALSE-----------------------------------------------------------------------
# ## selectVar(spls1.liver, comp = 1)$X$name
# 
# 
# ## ----------------------------------------------------------------------------------------
# spls1.liver$explained_variance$X
# pls1.liver.for.tune$explained_variance$X
# 
# 
# ## ----spls1-ext, fig.cap='(ref:spls1-ext)'------------------------------------------------
# spls1.liver.c2 <- spls(X, y, ncomp = 2, keepX = c(rep(choice.keepX, 2)), 
#                    mode = "regression")
# 
# plotIndiv(spls1.liver.c2,
#           group = liver.toxicity$treatment$Time.Group,
#           pch = as.factor(liver.toxicity$treatment$Dose.Group),
#           legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')
# 
# 
# 
# ## ----spls1-comp1, fig.cap='(ref:spls1-comp1)'--------------------------------------------
# plot(spls1.liver$variates$X[,1], spls1.liver$variates$Y[,1], 
#      xlab = 'X component', ylab = 'y component / scaled y',
#      main = 'Dimension 1')
# 
# cor(spls1.liver$variates$X, spls1.liver$variates$Y)
# 
# 
# ## ----------------------------------------------------------------------------------------
# set.seed(33)  # for reproducibility with this handbook, remove otherwise
# # PLS1 model and performance
# pls1.liver <- pls(X, y, ncomp = choice.ncomp, mode = "regression")
# perf.pls1.liver <- perf(pls1.liver, validation = "Mfold", folds =10, 
#                    nrepeat = 10, progressBar = FALSE)
# 
# 
# # sPLS1 performance
# perf.spls1.liver <- perf(spls1.liver, validation = "Mfold", folds = 10, 
#                    nrepeat = 10, progressBar = FALSE)
# 
# perf.pls1.liver$MSEP; perf.spls1.liver$MSEP
# 
# 
# ## ----------------------------------------------------------------------------------------
# dim(Y)
# 
# 
# ## ----pls2-Q2, fig.cap='(ref:pls2-Q2)'----------------------------------------------------
# pls2.liver.for.tune = pls(X = X, Y = Y, ncomp = 10, mode = 'regression')
# set.seed(33)  # for reproducibility with this handbook, remove otherwise
# Q2.pls2.liver = perf(pls2.liver.for.tune, validation = 'Mfold', folds = 10)
# 
# 
# plot(Q2.pls2.liver$Q2.total, xlab = 'Components', ylab = 'Q2 total')
# abline(h = 0.0975); abline(h = 0, col = 'red')

# save.image(file = 'buildignore/devel/handbook/env.RData')
load('buildignore/devel/handbook/env.RData')
X <- X[,1:100]
## ---- eval = FALSE-----------------------------------------------------------------------
list.keepX = c(seq(5, 50, 5), seq(60, 150, 10))
list.keepX = c(2, 5)
list.keepY = c(3:10)
list.keepY = c(3, 8)
ncomp = 2
nrepeat = 3


# update with tuning code later + outputs
set.seed(33)  # for reproducibility with this handbook, remove otherwise
sPLS.tune.reg.cor.liver = tune.spls(X, Y, test.keepX = list.keepX, test.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'regression', measure.tune = 'cor')
plot(sPLS.tune.reg.cor.liver)
# sPLS.tune.reg.cor.liver.ref <- sPLS.tune.reg.cor.liver
# PLS.tune.reg.cor.liver.ref <- PLS.tune.reg.cor.liver
# save(sPLS.tune.reg.cor.liver.ref, PLS.tune.reg.cor.liver.ref, file='res.RData')
load('res.RData')
# testthat::expect_equal(sPLS.tune.reg.cor.liver, sPLS.tune.reg.cor.liver.ref)
# testthat::expect_equal(PLS.tune.reg.cor.liver, PLS.tune.reg.cor.liver.ref)
waldo::compare(y = sPLS.tune.reg.cor.liver, sPLS.tune.reg.cor.liver.ref)
waldo::compare(y = PLS.tune.reg.cor.liver, PLS.tune.reg.cor.liver.ref)

testthat::expect_condition(tune.spls(X, Y, test.keepX = list.keepX, test.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'regression', measure.tune = 'cor', method = 'pls'))
sPLS.tune.reg.cor.liver$best.keepX
#  20 110
sPLS.tune.reg.cor.liver$best.keepY
# 3 4


## ----------------------------------------------------------------------------------------
choice.ncomp <- 2
choice.keepX <- c(20, 110)
choice.keepY <- c(3, 4)

liver.spls2 <- spls(X, Y, ncomp = choice.ncomp, 
                   keepX = choice.keepX,
                   keepY = choice.keepY,
                   mode = "regression")


## ----------------------------------------------------------------------------------------
liver.spls2$explained_variance


## ---- eval = FALSE-----------------------------------------------------------------------
## selectVar(liver.spls2, comp = 1)$X$value


## ----------------------------------------------------------------------------------------
liver.spls2.vip <- vip(liver.spls2)
liver.spls2.vip[selectVar(liver.spls2, comp = 1)$X$name,1]


## ----------------------------------------------------------------------------------------
liver.spls2.perf <- perf(liver.spls2, validation = 'Mfold', folds = 10)
liver.spls2.perf$features$stable.X$comp1[1:20]


## ----spls2-plotIndiv, fig.cap='(ref:spls2-plotIndiv)'------------------------------------
plotIndiv(liver.spls2, ind.names = FALSE, group = liver.toxicity$treatment$Time.Group, 
          pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')


## ----spls2-plotArrow, fig.cap='(ref:spls2-plotArrow)'------------------------------------
col.arrow = color.mixo(as.numeric(as.factor(liver.toxicity$treatment$Time.Group)))
plotArrow(liver.spls2, ind.names = FALSE, group = col.arrow)


## ----spls2-plotVar, fig.cap='(ref:spls2-plotVar)'----------------------------------------
plotVar(liver.spls2, cex = c(3,4), var.names = c(FALSE, TRUE))


## ----spls2-plotVar2, fig.cap='(ref:spls2-plotVar2)'--------------------------------------
plotVar(liver.spls2,
        var.names = list(X.label = liver.toxicity$gene.ID[,'geneBank'],
        Y.label = TRUE), cex = c(3,4))


## ----------------------------------------------------------------------------------------
liver.toxicity$gene.ID[selectVar(liver.spls2, comp = 1)$X$name,'geneBank']


## ----------------------------------------------------------------------------------------
# define red and green colors for the edges
color.edge <- color.GreenRed(50)
# X11()  # to open a new window if Rstudio complains
network(liver.spls2, comp = 1:2,
        cutoff = 0.7,
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        # to save the plot, otherwise comment out:
        save = 'pdf', name.save = 'network_liver')


## ----spls2-network, eval=TRUE, echo=FALSE,fig.cap='(ref:spls2-network)'------------------
knitr::include_graphics("network_liver.pdf")


## ----------------------------------------------------------------------------------------
# X11()  # to open a new window if Rstudio complains
cim(liver.spls2, comp = 1:2, xlab = "clinic", ylab = "genes",
    # to save the plot, otherwise comment out:
    save = 'pdf', name.save = 'cim_liver')


## ----spls2-cim, eval=TRUE, echo=FALSE, fig.cap='(ref:spls2-cim)'-------------------------
knitr::include_graphics("cim_liver.pdf")


## ---- eval = FALSE-----------------------------------------------------------------------
## # code to run and update here
## 
## # Comparisons of final models (PLS, sPLS)
## PLS.liver = pls(X, Y, mode = 'regression', ncomp = 3)
## source('perf_spls_repeat.R')
## perf.PLS.liver =  perf_spls_repeat(PLS.liver, validation = 'Mfold', folds = 10, nrepeat = 10)
## 
## sPLS.liver.cor = spls(X, Y, keepX = c(20,110), keepY = c(3,4), mode = 'regression', ncomp = 3)
## perf.sPLS.liver.cor =  perf_spls_repeat(sPLS.liver.cor, validation = 'Mfold', folds = 10, nrepeat = 10)
## 
## 
## plot(perf.PLS.liver$cor.pred$u, col = 'blue', pch = 16, ylim = c(0.1,1), ylab = 't or u Cor', main = 'PLS based on Cor')
## points(perf.PLS.liver$cor.pred$t, col = 'red', pch = 16)
## points(perf.sPLS.liver.cor$cor.pred$u, col = 'blue', pch = 17)
## points(perf.sPLS.liver.cor$cor.pred$t, col = 'red', pch = 17)
## legend('topright', col = c('blue', 'red', 'blue', 'red'), pch = c(16, 16, 17, 17), c('u PLS', 't PLS', 'u sPLS', 't sPLS'), cex = 0.5)


## ----pls-perf-spls2, eval=TRUE, echo=FALSE,  fig.cap='(ref:pls-perf-spls2)', out.width = '70%'----
knitr::include_graphics("Figures/Part3/liver_perf_spls2_reg.jpeg")


## ----------------------------------------------------------------------------------------
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
dim(X)
dim(Y)


## ----------------------------------------------------------------------------------------
linn.pls <- pls(X, Y, ncomp = 2, mode = "regression", scale = TRUE)


## ----------------------------------------------------------------------------------------
# new test individuals
indiv1 <- c(14, 190, 35)
indiv2 <- c(2, 40, 45)

X_test <- rbind(indiv1, indiv2)
colnames(X_test) <- colnames(X)


## ----------------------------------------------------------------------------------------
Y_pred <- predict(linn.pls, X_test)


## ----------------------------------------------------------------------------------------
Y_pred$predict


## ----------------------------------------------------------------------------------------
# first plot the individuals from the learning set in the X-space'
# graphics style needed here to overlay plots
plotIndiv(linn.pls, style = 'graphics', rep.space = 'X-variate',
          title =  'Linnerud PLS regression: sample plot')

# then the test set individuals coordinates that are predicted
points(Y_pred$variates[, 1], Y_pred$variates[, 2], pch = 17,
       col = c('red', 'blue'))
text(Y_pred$variates[, 1], Y_pred$variates[, 2], c("ind1", "ind2"), pos = 1,
     col = c('red', 'blue'))


## ---- eval = FALSE-----------------------------------------------------------------------
## data('nutrimouse')
## X = nutrimouse$gene
## Y = nutrimouse$lipid
## dim(X); dim(Y)
## 
## nrepeat = 10
## ncomp = 5
## list.keepX =  ncol(X)
## list.keepY = ncol(Y)
## 
## PLS.tune.can2 = tune_spls2_repeat(X, Y, list.keepX = list.keepX, list.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'canonical', type.tune = NULL, pls.model = TRUE)
## 
## PLS.tune.reg2 = tune_spls2_repeat(X, Y, list.keepX = list.keepX, list.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'regression', type.tune = NULL, pls.model = TRUE)
## 
## Q2.nutri.can2 = PLS.tune.can2$Q2.tot.ave
## Q2.nutri.reg2 = PLS.tune.reg2$Q2.tot.ave
## 
## 
## par(mfrow=c(1,2))
## plot(Q2.nutri.can2, xlab = 'Components', ylab = 'Q2 total Canonical', ylim = c(-0.2, 0.5))
## abline(h = 0.0975); abline(h = 0, col = 'red')
## title(main = 'Nutrimouse')
## 
## plot(Q2.nutri.reg2, xlab = 'Components', ylab = 'Q2 total Reg', ylim = c(-0.2, 0.5))
## abline(h = 0.0975); abline(h = 0, col = 'red')


## ----pls-canonical, eval=TRUE, echo=FALSE, fig.cap='(ref:pls-canonical)', fig.align="center", out.width = '70%'----
knitr::include_graphics("Figures/Part3/nutrimouse_Q2_PLS.pdf")


## ---- eval = FALSE-----------------------------------------------------------------------
## list.keepX = c(seq(5, 100, 5))
## list.keepY = c(seq(5, 20, 5))
## ncomp = 3
## nrepeat = 10
## 
## # may take a while to compute
## 
## source('tune_spls2_repeat.R')
## 
## sPLS.tune.can.cor.nutri = tune_spls2_repeat(X, Y, list.keepX = list.keepX, list.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'canonical', type.tune = 'cor', pls.model = FALSE)
## 
## sPLS.tune.can.RSS.nutri = tune_spls2_repeat(X, Y, list.keepX = list.keepX, list.keepY = list.keepY, ncomp = ncomp, nrepeat = nrepeat, mode = 'canonical', type.tune = 'RSS', pls.model = FALSE)
## 
## sPLS.tune.can.cor.nutri$best.keepX
## sPLS.tune.can.cor.nutri$best.keepY
## 
## sPLS.tune.can.RSS.nutri$best.keepX
## sPLS.tune.can.RSS.nutri$best.keepY


## ----pls-loc-reg, eval=TRUE, echo=FALSE, fig.cap='(ref:pls-loc-reg)', fig.align="center", out.width = '70%'----
knitr::include_graphics("Figures/Part3/PLS-loc-reg.pdf")

