## ----------- .get_xy block -----------
## example
## Y_mat ~ X1_mat + X2_mat
# mc = list(formula=Y_bpls ~ X_bpls$mirna + X_bpls$mrna)
# mf <- model.frame(mc$formula)
# mc$Y <- mf[1]
# mc$X <- as.list(mf[-1])
# mc$Y <- mc$formula[[2]]
# done

mc = as.call(list(block.pls, formula = protein ~ mirna + mrna, data = breast.TCGA.mae))

.get_xy_block <- function(mc) {
  trms <- rownames(attr(terms(mc$formula), "factors"))
  missing.assays <-  trms[!trms %in% names(experiments(mc$data))]
  if (length(missing.assays))
    .stop(paste("Not valid assay name(s) from data: ", paste(missing.assays, collapse = ", ")), .subclass = "inv_XY")
  ## get Y
  mc$Y <- assay(mc$data, trms[[1]])
  ## get X
  mc$X <- list()
  for (asy in trms[-1]) {
    mc$X[[asy]] <- assay(mc$data, asy)
  }
  mc$data <- mc$formula <- NULL
  mc
}