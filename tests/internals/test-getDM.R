library(testthat)
library(MultiAssayExperiment)

## with a valid assay name we get a matrix
expect_true(is.matrix(.getDM(X=miniACC, Assay = names(assays(miniACC))[1])))
expect_true(is.matrix(.getDM(X=miniACC, Assay = 1)))
## with an invalid assay name we get error
expect_error(.getDM(X=miniACC, Assay = NA), class = "invalid_assay")
expect_error(.getDM(X=miniACC, Assay = "invalidAssay"), class = "invalid_assay")
expect_error(.getDM(X=miniACC, Assay = 21), class = "invalid_assay")
expect_error(.getDM(X=miniACC, Assay = 20), class = "invalid_assay")
