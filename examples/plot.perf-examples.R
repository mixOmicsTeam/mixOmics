## ---- plot ffor perf.pls1 and perf.pls2
 data("nutrimouse")
 X <- data.matrix(nutrimouse$gene, rownames.force = TRUE)
 Y <- data.matrix(nutrimouse$lipid, rownames.force = TRUE)
 y <- Y[, "C18.0"]
 
 pls1.liver <- pls(X = X, Y = y, ncomp = 4)
 pls2.liver <- pls(X = X, Y = Y, ncomp = 4)
 spls1.liver <- spls(X = X, Y = y, keepX = rep(10, 4), ncomp = 4)
 spls2.liver <- spls(X = X, Y = Y, keepX = rep(10, 4), ncomp = 4)
 

 perf.pls1.liver  = perf(pls1.liver, validation = 'Mfold', 
                         folds = 3, nrepeat = 3, progressBar = TRUE)
 plot(perf.pls1.liver , criterion = 'Q2')
 perf.pls2.liver  = perf(pls2.liver, validation = 'Mfold', 
                         folds = 3, nrepeat = 3, progressBar = TRUE)
 plot(perf.pls2.liver , criterion = 'Q2.total')
 perf.spls1.liver = perf(spls1.liver, validation = 'Mfold', 
                         folds = 3, nrepeat = 3, progressBar = TRUE)
 plot(perf.spls1.liver , criterion = 'Q2')
 perf.spls2.liver = perf(spls2.liver, validation = 'Mfold', 
                         folds = 3, nrepeat = 3, progressBar = TRUE)
 plot(perf.spls2.liver , criterion = 'Q2.total')
 
 ## ---- plot for perf.(s)plsda
\dontrun{
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
set.seed(42)
breast.perf <- perf(res, nrepeat = 5)

plot(breast.perf)
plot(breast.perf, col=1:3)
plot(breast.perf, col=1:3, sd=FALSE)
}