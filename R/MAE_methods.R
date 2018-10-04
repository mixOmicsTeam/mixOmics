################################################################################
# Authors:
#   Sebastien Dejean
#   Florian Rohart
#
# created: 27-09-2018
# last modified: 27-09-2018
#
# Copyright (C) 2018
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
################################################################################

if(FALSE){
setGeneric(name = "pca", def = function(X, Y=NULL, ...) standardGeneric("pca"))

setMethod(f = "pca",
signature = signature("MultiAssayExperiment","character"),
definition = function(X, Y, ...) {
    tdm <- function(x) data.matrix(t(x))
    mixOmics::pca( tdm(experiments(X)[[Y]]), ... ) } )

setMethod(f = "pca",  # with a one-sided formula ~lipid
signature = signature("MultiAssayExperiment","formula"),
definition = function(X, Y, ...) {
    Y = as.character(as.list(Y)[[2]])
    tdm <- function(x) data.matrix(t(x))
    mixOmics::pca( tdm(experiments(X)[[Y]]), ... )    
} )

setMethod(f = "pca",
signature = "data.frame",
definition = function(X, Y=NULL, ...)  mixOmics::pca(X, ...) )


setGeneric(name = "spca", def = function(X, Y=NULL,  ...) standardGeneric("spca"))

setMethod(f = "spca",
signature = signature("MultiAssayExperiment", "character"),
definition = function(X, Y, ...) {
    tdm <- function(x) data.matrix(t(x))
    mixOmics::spca( tdm(experiments(X)[[Y]]), ...)
})


setMethod(f = "spca",
signature = signature("MultiAssayExperiment","formula"),
definition = function(X, Y, ...) {
    Y = as.character(as.list(Y)[[2]])
    tdm <- function(x) data.matrix(t(x))
    mixOmics::spca( tdm(experiments(X)[[Y]]), ... ) } )

setMethod(f = "spca", signature=signature("data.frame"),
definition = function(X, Y=NULL, ...) mixOmics::spca(X, ...))



setGeneric(name = "ipca", def = function(X, Y=NULL, ...) standardGeneric("ipca"))

setMethod(f = "ipca",
signature = signature("MultiAssayExperiment","character"),
definition = function(X, Y, ...) {
    tdm <- function(x) data.matrix(t(x))
    mixOmics::ipca( tdm(experiments(X)[[Y]]), ... ) } )

setMethod(f = "ipca",
signature = signature("MultiAssayExperiment","formula"),
definition = function(X, Y, ...) {
    Y = as.character(as.list(Y)[[2]])
    tdm <- function(x) data.matrix(t(x))
    mixOmics::ipca( tdm(experiments(X)[[Y]]), ... ) } )

setMethod(f = "ipca",
signature = "data.frame",
definition = function(X, Y=NULL, ...)  mixOmics::ipca(X, ...) )



setGeneric(name = "sipca", def = function(X, Y=NULL,  ...) standardGeneric("sipca"))

setMethod(f = "sipca",
signature = signature("MultiAssayExperiment", "character"),
definition = function(X, Y, ...) {
    tdm <- function(x) data.matrix(t(x))
    mixOmics::sipca( tdm(experiments(X)[[Y]]), ...)
})

setMethod(f = "sipca",
signature = signature("MultiAssayExperiment","formula"),
definition = function(X, Y, ...) {
    Y = as.character(as.list(Y)[[2]])
    tdm <- function(x) data.matrix(t(x))
    mixOmics::sipca( tdm(experiments(X)[[Y]]), ... ) } )

setMethod(f = "sipca", signature=signature("data.frame"),
definition = function(X, Y=NULL, ...) mixOmics::sipca(X, ...))



setGeneric(name = "plsda", def = function(X, Y, ...) standardGeneric("plsda"))

setMethod(f = "plsda",
signature = signature("MultiAssayExperiment", "formula"),
definition = function(X, Y, ...) {
    els <- vapply(as.list(Y), as.character, "character")[-1]
    tdm <- function(x) data.matrix(t(x))
    ans <- mixOmics::plsda(tdm(experiments(X)[[els[2]]]),
    colData(X)[[els[1]]], ...) } )

setMethod(f = "plsda", signature=signature("data.frame", "factor"),
definition = function(X, Y, ...) mixOmics::plsda(X, Y, ...))



setGeneric(name = "splsda", def = function(X, Y, ...) standardGeneric("splsda"))

setMethod(f = "splsda",
signature = signature("MultiAssayExperiment", "formula"),
definition = function(X, Y, ...) {
    els <- vapply(as.list(Y), as.character, "character")[-1]
    tdm <- function(x) data.matrix(t(x))
    ans <- mixOmics::splsda(tdm(experiments(X)[[els[2]]]),
    colData(X)[[els[1]]], ...) } )

setMethod(f = "splsda", signature=signature("data.frame", "factor"),
definition = function(X, Y, ...) mixOmics::splsda(X, Y, ...))



setGeneric("pls", function(X, Y, ...) standardGeneric("pls"))

# MAE + formula => mode regression
setMethod("pls", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula: regression mode ***")
    els <- vapply(as.list(Y), as.character, "character")[-1]
    tdm <- function(x) data.matrix(t(x))
    mixOmics::pls( tdm(experiments(X)[[els[1]]]),
    tdm(experiments(X)[[els[2]]]), mode="regression", ...)
})

# MAE + vector(character + character) => mode canonical
setMethod("pls", c("MultiAssayExperiment", "vector"),
function(X, Y, ...) {
    message("*** MAE + 2-vector of characters ***")
    tdm <- function(x) data.matrix(t(x))
    mixOmics::pls( tdm(experiments(X)[[Y[1]]]),
    tdm(experiments(X)[[Y[2]]]), ...)
})

# A l'ancienne : 2 data.frame
setMethod("pls", c("data.frame", "data.frame"), function(X, Y, ...) {
    message("*** Two data.frames ***")
    mixOmics::pls(X, Y, ...)
})




setGeneric("spls", function(X, Y, ...) standardGeneric("spls"))

# MAE + formula => mode regression
setMethod("spls", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula: regression mode ***")
    els <- vapply(as.list(Y), as.character, "character")[-1]
    tdm <- function(x) data.matrix(t(x))
    mixOmics::spls( tdm(experiments(X)[[els[1]]]),
    tdm(experiments(X)[[els[2]]]), mode="regression", ...)
})

# MAE + vector(character + character) => mode canonical
setMethod("spls", c("MultiAssayExperiment", "vector"),
function(X, Y, ...) {
    message("*** MAE + 2-vector of characters ***")
    tdm <- function(x) data.matrix(t(x))
    mixOmics::spls( tdm(experiments(X)[[Y[1]]]),
    tdm(experiments(X)[[Y[2]]]), ...)
})

# A l'ancienne : 2 data.frame
setMethod("spls", c("data.frame", "data.frame"), function(X, Y, ...) {
    message("*** Two data.frames ***")
    mixOmics::spls(X, Y, ...)
})


setGeneric("rcc", function(X, Y, ...) standardGeneric("rcc"))

# MAE + vector(character + character) => mode canonical
setMethod("rcc", c("MultiAssayExperiment", "vector"),
function(X, Y, ...) {
    message("*** MAE + 2-vector of characters ***")
    tdm <- function(x) data.matrix(t(x))
    mixOmics::rcc( tdm(experiments(X)[[Y[1]]]),
    tdm(experiments(X)[[Y[2]]]), ...)
})

# A l'ancienne : 2 data.frame
setMethod("rcc", c("data.frame", "data.frame"), function(X, Y, ...) {
    message("*** Two data.frames ***")
    mixOmics::rcc(X, Y, ...)
})

# No formula with MAE because it could be misleading



setGeneric("block.pls", function(X, Y, ...) standardGeneric("block.pls"))

setMethod("block.pls", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula ***")
    tdm <- function(x) data.matrix(t(x))
    list_sets <- names(experiments(X))
    els <- as.character(Y)[-1]
    Y.input <- tdm(experiments(X)[[els[1]]])
    if (els[2] == ".")
    {
        list_X_input <- setdiff(list_sets, els[1])
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    else
    {
        list_X_input <- unlist(strsplit(els[2],"[ + ]"))
        list_X_input <- list_X_input[list_X_input != ""]
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    mixOmics::block.pls(X.input, Y.input , ...)
})

# A l'ancienne : X list, Y data.frame
setMethod("block.pls", c("list", "data.frame"),
function(X, Y, ...) {
    message("*** One list + one data.frame ***")
    mixOmics::block.pls(X, Y, ...)
})

# ou X list, Y matrix
setMethod("block.pls", c("list", "matrix"),
function(X, Y, ...) {
    message("*** One list + one matrix ***")
    mixOmics::block.pls(X, Y, ...)
})



setGeneric("block.spls", function(X, Y, ...) standardGeneric("block.spls"))

setMethod("block.spls", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula ***")
    tdm <- function(x) data.matrix(t(x))
    list_sets <- names(experiments(X))
    els <- as.character(Y)[-1]
    Y.input <- tdm(experiments(X)[[els[1]]])
    if (els[2] == ".")
    {
        list_X_input <- setdiff(list_sets, els[1])
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    else
    {
        list_X_input <- unlist(strsplit(els[2],"[ + ]"))
        list_X_input <- list_X_input[list_X_input != ""]
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    mixOmics::block.spls(X.input, Y.input , ...)
})

# A l'ancienne : X list, Y data.frame
setMethod("block.spls", c("list", "data.frame"),
function(X, Y, ...) {
    message("*** One list + one data.frame ***")
    mixOmics::block.spls(X, Y, ...)
})

# ou X list, Y matrix
setMethod("block.spls", c("list", "matrix"),
function(X, Y, ...) {
    message("*** One list + one matrix ***")
    mixOmics::block.spls(X, Y, ...)
})



setGeneric("block.plsda", function(X, Y, ...) standardGeneric("block.plsda"))

setMethod("block.plsda", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula ***")
    tdm <- function(x) data.matrix(t(x))
    list_sets <- names(experiments(X))
    els <- as.character(Y)[-1]
    Y.input <- as.factor(colData(X)[els[1]][,1])
    if (els[2] == ".")
    {
        list_X_input <- list_sets
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    else
    {
        list_X_input <- unlist(strsplit(els[2],"[ + ]"))
        list_X_input <- list_X_input[list_X_input != ""]
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    mixOmics::block.plsda(X.input, Y.input , ...)
})

# A l'ancienne : X list, Y factor
setMethod("block.plsda", c("list", "factor"),
function(X, Y, ...) {
    message("*** One list + one factor ***")
    mixOmics::block.plsda(X, Y, ...)
})



setGeneric("block.splsda", function(X, Y, ...) standardGeneric("block.splsda"))

setMethod("block.splsda", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula ***")
    tdm <- function(x) data.matrix(t(x))
    list_sets <- names(experiments(X))
    els <- as.character(Y)[-1]
    Y.input <- as.factor(colData(X)[els[1]][,1])
    if (els[2] == ".")
    {
        list_X_input <- list_sets
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    else
    {
        list_X_input <- unlist(strsplit(els[2],"[ + ]"))
        list_X_input <- list_X_input[list_X_input != ""]
        ind <- length(list_X_input)
        X.input <- list()
        for (i in 1:ind) {
            X.input[[i]] <- tdm(experiments(X)[[list_X_input[i]]])
            names(X.input)[i] <- list_X_input[i]
        }
    }
    mixOmics::block.splsda(X.input, Y.input , ...)
})

# A l'ancienne : X list, Y factor
setMethod("block.splsda", c("list", "factor"),
function(X, Y, ...) {
    message("*** One list + one factor ***")
    mixOmics::block.splsda(X, Y, ...)
})



setGeneric("wrapper.rgcca", function(X, Y = NULL, ...) standardGeneric("wrapper.rgcca"))

# Si rien n'est pr?cis? en Y, on garde tous les experiments du MAE
setMethod("wrapper.rgcca", "MultiAssayExperiment",
function(X, Y = NULL, ...) {
    message("*** MAE only ***")
    tdm <- function(x) data.matrix(t(x))
    list_sets <- names(experiments(X))
    ind <- length(list_sets)
    X.input <- list()
    for (i in 1:ind) {
        X.input[[i]] <- tdm(experiments(X)[[list_sets[i]]])
        names(X.input)[i] <- list_sets[i]
    }
    mixOmics::wrapper.rgcca(X.input, ...)
})

# On peut donner en Y les experiments ? prendre en compte dans un vecteur de character
setMethod("wrapper.rgcca", c("MultiAssayExperiment", "vector"),
function(X, Y, ...) {
    message("*** MAE + vector ***")
    tdm <- function(x) data.matrix(t(x))
    ind <- length(Y)
    X.input <- list()
    for (i in 1:ind) {
        X.input[[i]] <- tdm(experiments(X)[[Y[i]]])
        names(X.input)[i] <- Y[i]
    }
    mixOmics::wrapper.rgcca(X.input, ...)
})

# A l'ancienne : X list
setMethod("wrapper.rgcca", "list",
function(X, Y = NULL, ...) {
    message("*** One list ***")
    mixOmics::wrapper.rgcca(X, ...)
})


setGeneric("wrapper.sgcca", function(X, Y = NULL, ...) standardGeneric("wrapper.sgcca"))

# Si rien n'est pr?cis? en Y, on garde tous les experiments du MAE
setMethod("wrapper.sgcca", "MultiAssayExperiment",
function(X, Y = NULL, ...) {
    message("*** MAE only ***")
    tdm <- function(x) data.matrix(t(x))
    list_sets <- names(experiments(X))
    ind <- length(list_sets)
    X.input <- list()
    for (i in 1:ind) {
        X.input[[i]] <- tdm(experiments(X)[[list_sets[i]]])
        names(X.input)[i] <- list_sets[i]
    }
    mixOmics::wrapper.sgcca(X.input, ...)
})

# On peut donner en Y les experiments ? prendre en compte dans un vecteur de character
setMethod("wrapper.sgcca", c("MultiAssayExperiment", "vector"),
function(X, Y, ...) {
    message("*** MAE + vector ***")
    tdm <- function(x) data.matrix(t(x))
    ind <- length(Y)
    X.input <- list()
    for (i in 1:ind) {
        X.input[[i]] <- tdm(experiments(X)[[Y[i]]])
        names(X.input)[i] <- Y[i]
    }
    mixOmics::wrapper.sgcca(X.input, ...)
})

# A l'ancienne : X list
setMethod("wrapper.sgcca", "list",
function(X, Y = NULL, ...) {
    message("*** One list ***")
    mixOmics::wrapper.sgcca(X, ...)
})



setGeneric("wrapper.sgccda", function(X, Y, ...) standardGeneric("wrapper.sgccda"))

setMethod("wrapper.sgccda", c("MultiAssayExperiment", "formula"),
function(X, Y, ...) {
    message("*** MAE + formula ***")
    block.splsda(X, Y , ...)
})

# Il faudra certainement adapter le truc en l'int?grant au package.
# Ici, je ne fais que reprendre la fonction block.plsda definie plus haut
# sans passer par la fonction block.plsda de l'actuel package mixOmics

# A l'ancienne : X list, Y factor
setMethod("wrapper.sgccda", c("list", "factor"),
function(X, Y, ...) {
    message("*** One list + one factor ***")
    block.splsda(X, Y, ...)
})

}
