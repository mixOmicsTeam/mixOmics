
## example with PLSDA, 2 classes
# ----------------
#' \dontrun{
data(breast.tumors)

X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
auc.plsda.breast = auroc(plsda.breast, ncomp = 1)


  ## example with sPLSDA
  # -----------------
  splsda.breast <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
  auroc(plsda.breast, plot = FALSE)


  ## example with sPLSDA with 4 classes
  # -----------------
  X <- as.matrix(liver.toxicity$gene)
  # Y will be transformed as a factor in the function,
  # but we set it as a factor to set up the colors.
  Y <- as.factor(liver.toxicity$treatment[, 4])

  splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))
  auc.splsda.liver = auroc(splsda.liver, ncomp = 1)


  ## example with mint.plsda
  # -----------------

  res = mint.plsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3,
                   study = stemcells$study)
  auc.mint.pslda = auroc(res, plot = FALSE)

  ## example with mint.splsda
  # -----------------
  data(stemcells)
  res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
                    study = stemcells$study)
  auc.mint.spslda = auroc(res, plot = TRUE, roc.comp = 3)


  ## example with block.plsda
  # ------------------
  data(nutrimouse)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
  # with this design, all blocks are connected
  design = matrix(c(0,1,1,0), ncol = 2, nrow = 2,
                  byrow = TRUE, dimnames = list(names(data), names(data)))

  block.plsda.nutri = block.plsda(X = data, Y = nutrimouse$diet)
  auc.block.plsda.nutri = auroc(block.plsda.nutri, block = 'lipid')

  ## example with block.splsda
  # ---------------
  list.keepX = list(gene = rep(10, 2), lipid = rep(5,2))
  block.splsda.nutri = block.splsda(X = data, Y = nutrimouse$diet, keepX = list.keepX)
  auc.block.splsda.nutri = auroc(block.splsda.nutri, block = 1)
#' }
