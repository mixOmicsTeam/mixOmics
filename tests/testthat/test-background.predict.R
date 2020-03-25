context("background.predict")

test_that("background.predict works", code = {
    data(liver.toxicity)
    X = liver.toxicity$gene
    Y = as.factor(liver.toxicity$treatment[, 4])
    
    plsda.liver <- plsda(X, Y, ncomp = 2)
    background = background.predict(plsda.liver, comp.predicted = 2, dist = "mahalanobis.dist", resolution = 20)
    
    expect_is(background, "background.predict")
    .expect_numerically_close(background$`6`[1,], c(Var1 = 16.3966070584067, Var2 = -28.7410902930419))
})
