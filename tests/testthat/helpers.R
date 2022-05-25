#' Expect identical ignoring $call slot
#'
#' @param object an R object
#' @param expected an R object
#' @param ... other args passed to expect_identical
#' @keywords internal
#' @return Logical, TRUE if all elements expect $call are identical
.almost_identical <- function(object, expected, ...) {
    object$call <- expected$call <- NULL
    expect_identical(object, expected, ...)
}


#' Expect equal numeric after rounding
#'
#' @param numeric_value numeric, outcome of test
#' @param expected numeric, reference to check
#' @param digits integer, decimals to round
#' @return logical, if TRUE, values are almost equal
#' @keywords internal
#' @examples
#' .numerically_close(3.414, 3.14, digits =2)
.expect_numerically_close <- function(numeric_value, expected, digits = 2) {
    require(testthat)
    expect_equal(round(numeric_value, digits = digits), round(expected, digits = digits))
}


#' From input X and Y dataframes, yields the smallest set of training and testing
#' samples to remain valid for any mixOmics method. Caters sample selection to
#' if method requires multiblock, multigroup or multilevel frameworks.
#' 
#' @param X X dataframe for any mixOmics method. Can be a list of multiple dataframes if multiblock
#' @param Y Y dataframe or factor vector for any mixOmics method
#' @param S study factor vector for multigroup frameworks
#' @param ML repreated measures vector for multilevel frameworks
#' @param seed controls the sample selection seed
#' @return list of X, Y, study and multilevel components split by training and testing samples
#' @keywords internal
.minimal_train_test_subset <- function(X=NULL, Y=NULL, S=NULL, ML=NULL,
                                       seed=16) {
    set.seed(seed)
    
    DA = is.factor(Y) # logical gate for DA framework
    MULTIGROUP = !is.null(S) # logical gate for multigroup framework
    MULTILEVEL = !is.null(ML) # logical gate for multilevel framework
    MULTIBLOCK = !is.data.frame(X) && !is.matrix(X) # logical gate for multiblock framework
    
    tr <- c() # initialise indicies for training and testing samples
    te <- c()
    
    if (MULTILEVEL) { # any multilevel method
        
        n.indivs <- 3 # default case for number of repeated samples to consider
        
        if(DA) { n.indivs <- length(unique(Y))-1 } # if DA, set specific quantity
            
        # only look at the first n.indiv samples were measured the maximum amount of times
        indivs <- unname(which(table(ML) == max(table(ML))))[1:n.indivs] 
        
        for (i in 1:length(indivs)) { # for each repeated sample ...
            s <- indivs[i]
            
            rel.sam <- which(ML==s) # determine the corresponding rows
            tr.sam <- rel.sam[1:2+(i-1)] # take two of these for training
            te.sam <- setdiff(rel.sam, tr.sam) # and the remaining for testing
            
            te <- c(te, tr.sam)
            tr <- c(tr, te.sam)
            
        }
    } 
    else if(DA) { # if the framework is DA ... 
        
        for(c in unique(Y)) { # for each class ...
            
            if (MULTIGROUP) { # MINT.(s)PLSDA
                for (s in unique(S)){ # for each study ...
                    # determine the rows with that class and for that study
                    rel.sam <- intersect(which(Y==c), which(S==s)) 
                    tr <- c(tr, rel.sam[1]) # take first for training
                    # if that samples's class and study is not already present in testing, add it
                    if (!(s %in% S[te] || c %in% Y[te])) {te <- c(te, rel.sam[2]) }
                }
            } else { # (BLOCK).(s)PLSDA
                rows <- which(Y == c)
                tr <- c(tr, rows[2:3])
                te <- c(te, rows[1])
            }
            
        }
        
        if (MULTIGROUP) { # ensure that all studies in training are present in testing
            tr.te.study.diff <- setdiff(unique(S[tr]), unique(S[te]))
            if (length(tr.te.study.diff) != 0) {
                for (s in tr.te.study.diff) {
                    te <- c(te, which(S == s)[1])
                }
            }
        }
        
    } 
    else { 
        if (MULTIGROUP) { # MINT.(S)PLS
            for (s in unique(S)){
                rel.sam <- which(S==s)
                te <- c(te, rel.sam[1])
                tr <- c(tr, rel.sam[c(2,3)])
            }
        } else { # (BLOCK).(s)PLS
            tr <- 1:6
            te <- 7:9
        }
    }
    
    
    
    if(MULTIBLOCK) { # subset each block iteratively if multiblock
        X.tr <- list()         
        X.te <- list()
        
        for (block in names(X)) {
            X.tr[[block]] <- X[[block]][tr,]
            X.te[[block]] <- X[[block]][te,]
        }
    } else { # otherwise just subset X
        X.tr <- X[tr, ]
        X.te <- X[te, ]
    }
    
    if (DA) { # if Y is a factor, index list
        Y.tr <- Y[tr]
        Y.te <- Y[te]
    } else { # if Y is a data.frame, index df
        Y.tr <- Y[tr,]
        Y.te <- Y[te,]
    }
    
    out <- list(X.tr = X.tr,
                X.te = X.te,
                Y.tr = Y.tr,
                Y.te = Y.te)
    
    if (MULTILEVEL) { # include repeated measures
        out$ML.tr <- ML[tr]
        out$ML.te <- ML[te]
    }
    
    if (MULTIGROUP) { # include study
        out$S.tr <- S[tr]
        out$S.te <- S[te]
    }
    
    
    return(out)
}
