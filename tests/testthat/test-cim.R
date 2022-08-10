
###############################################################################
### ============================ MISSING TESTS ============================ ###
###############################################################################

# basic
## pca - once issue #171 is resolved we can do this one

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################


test_that("(cim:basic): matrices", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid[, 1:10]
    Y <- nutrimouse$gene[, 1:10]
    
    matrix.cim <- cim(cor(X, Y), cluster = "none")
    
    expect_equal(class(matrix.cim), "cim_default")
    expect_equal(matrix.cim$row.names[4], "C16.1n.9")
    expect_equal(matrix.cim$col.names[4], "ACBP")
    .expect_numerically_close(matrix.cim$mat[4,4],
                              -0.197)
})


test_that("(cim:basic): spca", {

    data(nutrimouse)
    X <- nutrimouse$lipid[1:10,1:10]

    res.spca <- spca(X, keepX = c(3,3))

    spca.cim <- cim(res.spca)
    
    expect_equal(class(spca.cim), "cim_spca")
    expect_equal(spca.cim$row.names[4], "2")
    expect_equal(spca.cim$col.names[4], "C14.0")
    .expect_numerically_close(spca.cim$mat[4,4],
                              -0.356)
})


test_that("(cim:basic): ipca", {

    data(nutrimouse)
    X <- nutrimouse$lipid[1:10,1:10]

    res.ipca <- ipca(X)

    ipca.cim <- cim(res.ipca)
    
    expect_equal(class(ipca.cim), "cim_ipca")
    expect_equal(ipca.cim$row.names[4], "2")
    .expect_numerically_close(ipca.cim$mat[4,4],
                              -1.851)
})


test_that("(cim:basic): sipca", {

    data(nutrimouse)
    X <- nutrimouse$lipid[1:10,1:10]

    res.sipca <- sipca(X, ncomp=2, keepX = c(3,3))

    sipca.cim <- cim(res.sipca)
    
    expect_equal(class(sipca.cim), "cim_sipca")
    expect_equal(sipca.cim$row.names[4], "9")
    expect_equal(sipca.cim$col.names[4], "C16.1n.7")
    .expect_numerically_close(sipca.cim$mat[4,4],
                              -1.386)
})


test_that("(cim:basic): rcc", {

    data(nutrimouse)
    X <- nutrimouse$lipid[,1:10]
    Y <- nutrimouse$gene[,1:10]

    res.rcc <- rcc(X, Y, method = "shrinkage")

    rcc.cim <- cim(res.rcc)
    
    expect_equal(class(rcc.cim), "cim_rcc")
    expect_equal(rcc.cim$row.names[4], "C20.1n.9")
    expect_equal(rcc.cim$col.names[4], "X36b4")
    .expect_numerically_close(rcc.cim$mat[4,4],
                              -0.065)
})


test_that("(cim:basic): pls", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.pls <- pls(X, Y)

    pls.cim <- cim(res.pls)
    
    expect_equal(class(pls.cim), "cim_mixo_pls")
    expect_equal(pls.cim$row.names[4], "A_43_P12356")
    expect_equal(pls.cim$col.names[4], "TBA.umol.L.")
    .expect_numerically_close(pls.cim$mat[4,4],
                              0.109)
})


test_that("(cim:basic): spls", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3))

    spls.cim <- cim(res.spls)
    
    expect_equal(class(spls.cim), "cim_mixo_spls")
    expect_equal(spls.cim$row.names[4], "A_43_P12356")
    expect_equal(spls.cim$col.names[4], "SDH.IU.L.")
    .expect_numerically_close(spls.cim$mat[4,4],
                              -0.114)
})


test_that("(cim:basic): plsda", {
    
    set.seed(16)
    samples <- sample(1:63, 10)

    data(srbct)
    X <- srbct$gene[samples,1:10]
    Y <- srbct$class[samples]

    res.plsda <- plsda(X, Y)

    plsda.cim <- cim(res.plsda)
    
    expect_equal(class(plsda.cim), "cim_mixo_plsda")
    expect_equal(plsda.cim$row.names[4], "NB.C8")
    expect_equal(plsda.cim$col.names[4], "g3")
    .expect_numerically_close(plsda.cim$mat[4,4],
                              -0.351)
})


test_that("(cim:basic): splsda", {

    set.seed(16)
    samples <- sample(1:63, 10)

    data(srbct)
    X <- srbct$gene[samples,1:10]
    Y <- srbct$class[samples]

    res.splsda <- splsda(X, Y, keepX = c(3,3))

    splsda.cim <- cim(res.splsda)
    
    expect_equal(class(splsda.cim), "cim_mixo_splsda")
    expect_equal(splsda.cim$row.names[4], "NB.C2")
    expect_equal(splsda.cim$col.names[4], "g6")
    .expect_numerically_close(splsda.cim$mat[4,4],
                              -0.526)
})


test_that("(cim:basic): ml-spls", {

    data(vac18)
    X <- vac18$genes[,1:10]
    Y <- vac18$genes[,11:20]
    ml <- vac18$sample

    res.mlpls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3), multilevel = ml)

    mlspls.cim <- cim(res.mlpls)
    
    expect_equal(class(mlspls.cim), "cim_mixo_mlspls")
    expect_equal(mlspls.cim$row.names[4], "ILMN_1658145")
    expect_equal(mlspls.cim$col.names[4], "ILMN_1724774")
    .expect_numerically_close(mlspls.cim$mat[4,4],
                              -0.411)
})


test_that("(cim:basic): ml-splsda", {
    
    data(vac18)
    X <- vac18$genes[,1:10]
    Y <- vac18$stimulation
    ml <- vac18$sample

    res.mlplsda <- splsda(X, Y, keepX = c(3,3), multilevel = ml)

    mlsplsda.cim <- cim(res.mlplsda)
    
    expect_equal(class(mlsplsda.cim), "cim_mixo_mlsplsda")
    expect_equal(mlsplsda.cim$row.names[4], "X.55")
    expect_equal(mlsplsda.cim$col.names[4], "ILMN_1800889")
    .expect_numerically_close(mlsplsda.cim$mat[4,4],
                              1.276)
})


test_that("(cim:basic): mint.spls", {

    data(stemcells)
    X <- stemcells$gene[, 1:10]
    Y <- stemcells$gene[, 1:10]
    s <- stemcells$study

    res.mint.spls <- mint.spls(X, Y, study=s, keepX = c(3,3), keepY = c(3,3))

    mint.spls.cim <- cim(res.mint.spls)
    
    expect_equal(class(mint.spls.cim), "cim_mint.spls")
    expect_equal(mint.spls.cim$row.names[4], "ENSG00000189367")
    expect_equal(mint.spls.cim$col.names[4], "ENSG00000189367")
    .expect_numerically_close(mint.spls.cim$mat[4,4],
                              0.629)
})


test_that("(cim:basic): mint.splsda", {
    
    set.seed(17)
    samples <- sample(1:125, 40)

    data(stemcells)
    X <- stemcells$gene[samples, 1:10]
    Y <- stemcells$celltype[samples]
    s <- c(rep(1,10), rep(2,10), rep(3,10), rep(4,10))

    res.mint.splsda <- mint.splsda(X, Y, study=s, keepX = c(3,3))

    mint.splsda.cim <- cim(res.mint.splsda)
    
    expect_equal(class(mint.splsda.cim), "cim_mint.splsda")
    expect_equal(mint.splsda.cim$row.names[4], "sample146")
    expect_equal(mint.splsda.cim$col.names[4], "ENSG00000159199")
    .expect_numerically_close(mint.splsda.cim$mat[4,4],
                              0.337)
})


# ###############################################################################
# ### ================================ DATA ================================= ###
# ###############################################################################


test_that("(cim:data): multidrug", {

    data(multidrug)
    X <- multidrug$ABC.trans[1:10, 1:10]

    res.multidrug <- spca(X, keepX = c(3,3))

    multidrug.cim <- cim(res.multidrug)
    
    expect_equal(class(multidrug.cim), "cim_spca")
    expect_equal(multidrug.cim$row.names[4], "5")
    expect_equal(multidrug.cim$col.names[4], "ABCA4")
    .expect_numerically_close(multidrug.cim$mat[4,4],
                              -0.253)
})


# ###############################################################################
# ### ============================= PARAMETER =============================== ###
# ###############################################################################


test_that("(cim:paramter): cutoff", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    cutoff.cim <- cim(res.spls, cutoff = 0.2)

    expect_equal(class(cutoff.cim), "cim_mixo_spls")
    expect_equal(cutoff.cim$row.names[2], "A_43_P15834")
    expect_equal(cutoff.cim$col.names[3], "AST.IU.L.")
    .expect_numerically_close(cutoff.cim$mat[2,3],
                              -0.212)
})


test_that("(cim:paramter): cut.tree", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    cut.tree.cim <- cim(res.spls, cut.tree = c(0.5, 0.5))

    expect_equal(class(cut.tree.cim), "cim_mixo_spls")
    expect_equal(cut.tree.cim$row.names[4], "A_43_P22018")
    expect_equal(cut.tree.cim$col.names[3], "AST.IU.L.")
    .expect_numerically_close(cut.tree.cim$mat[4,3],
                              0.071)
})


test_that("(cim:paramter): cluster", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    # cluster = none
    cluster.none.cim <- cim(res.spls, cluster = "none")
    .expect_numerically_close(cluster.none.cim$mat[4,3],
                              0.126)
    
    # cluster = row
    cluster.row.cim <- cim(res.spls, cluster = "row")
    .expect_numerically_close(cluster.row.cim$mat[4,3],
                              -0.00948)

    # cluster = col
    cluster.col.cim <- cim(res.spls, cluster = "col")
    .expect_numerically_close(cluster.col.cim$mat[4,3],
                              0.199)
})


test_that("(cim:paramter): dist.method", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    dist.correlation.cim <- cim(res.spls, dist.method = c("correlation", "correlation"),
                                cluster = "row")

    expect_equal(class(dist.correlation.cim), "cim_mixo_spls")
    expect_equal(dist.correlation.cim$row.names[4], "A_43_P22018")
    expect_equal(dist.correlation.cim$col.names[3], "TBA.umol.L.")
    .expect_numerically_close(dist.correlation.cim$mat[4,3],
                              -0.00948)

    dist.manhattan.cim <- cim(res.spls, dist.method = c("manhattan", "manhattan"),
                              cluster = "col")

    expect_equal(class(dist.manhattan.cim), "cim_mixo_spls")
    expect_equal(dist.manhattan.cim$row.names[4], "A_43_P12356")
    expect_equal(dist.manhattan.cim$col.names[3], "AST.IU.L.")
    .expect_numerically_close(dist.manhattan.cim$mat[4,3],
                              0.1993)
})


test_that("(cim:paramter): clust.method", {

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    clust.centroid.cim <- cim(res.spls, clust.method = c("centroid", "centroid"),
                                cluster = "row")
    clust.complete.cim <- cim(res.spls, clust.method = c("complete", "complete"),
                                cluster = "row")

    if (!all(clust.complete.cim$mat == clust.centroid.cim$mat)) { stop("Objects should be equal", call.=F) }

    expect_equal(class(clust.complete.cim), "cim_mixo_spls")
    expect_equal(clust.complete.cim$row.names[4], "A_43_P22018")
    expect_equal(clust.complete.cim$col.names[3], "TBA.umol.L.")
    .expect_numerically_close(clust.complete.cim$mat[4,3],
                              -0.00948)
})


test_that("(cim:paramter): comp",{

    data(liver.toxicity)
    X <- liver.toxicity$gene[,1:10]
    Y <- liver.toxicity$clinic[,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3), ncomp = 4)

    comp.cim <- cim(res.spls, comp = c(2,3))

    expect_equal(class(comp.cim), "cim_mixo_spls")
    expect_equal(comp.cim$row.names[4], "A_43_P22290")
    expect_equal(comp.cim$col.names[4], "SDH.IU.L.")
    .expect_numerically_close(comp.cim$mat[4,4],
                              0.00745)
})


test_that("(cim:paramter): mapping",{

    data(liver.toxicity)
    X <- liver.toxicity$gene[1:10,1:10]
    Y <- liver.toxicity$clinic[1:10,1:10]

    res.spls <- spls(X, Y, keepX = c(3,3), keepY = c(3,3))

    mapping.X.cim <- cim(res.spls, mapping = "X")

    expect_equal(class(mapping.X.cim), "cim_mixo_spls")
    expect_equal(mapping.X.cim$row.names[4], "ID211")
    expect_equal(mapping.X.cim$col.names[4], "A_43_P21286")
    .expect_numerically_close(mapping.X.cim$mat[4,4],
                              -0.713)


    mapping.Y.cim <- cim(res.spls, mapping = "Y")

    expect_equal(class(mapping.Y.cim), "cim_mixo_spls")
    expect_equal(mapping.Y.cim$row.names[4], "ID209")
    expect_equal(mapping.Y.cim$col.names[3], "BUN.mg.dL.")
    .expect_numerically_close(mapping.Y.cim$mat[4,3],
                              0.556)
})


###############################################################################
### ================================ ERROR ================================ ###
###############################################################################


test_that("(cim:error): cluster is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, cluster = "invalid.cluster"),
                 "cluster")
})


test_that("(cim:error): clust.method is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, clust.method = c("ward", "invalid.cluster")),
                 "clustering")
})


test_that("(cim:error): dist.method is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, dist.method = c("euclidean", "invalid.cluster")),
                 "distance")
})


test_that("(cim:error): cut.tree is valid value", {
    
    data(nutrimouse)
    X <- nutrimouse$lipid
    
    res.spca <- spca(X)
    
    expect_error(cim(res.spca, cut.tree=0.8),
                 "cut.tree")
    
    expect_error(cim(res.spca, cut.tree=c(2,2)),
                 "cut.tree")
})


test_that("(cim:error): doesn't allow block.splsda objects", {
    
    data(breast.TCGA)
    X = list(miRNA = breast.TCGA$data.train$mirna,
             mRNA = breast.TCGA$data.train$mrna,
             proteomics = breast.TCGA$data.train$protein)
    Y = breast.TCGA$data.train$subtype
    
    res.block.splsda <- block.splsda(X, Y)
    
    expect_error(cim(res.block.splsda),
                 "cimDiablo")
})