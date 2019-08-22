## ----------- .plsMethodsHelper ----------- 
####  get call list including potentiall X, Y, formula, and data and retain only valid X and Y
.plsMethodsHelper <- function(mc){
  mc[c('data', 'formula')] <- lapply( mc[c('data', 'formula')], eval.parent)
  mcc <- mc ## copy so can change mc but keep the call for .get_xy
  # expectedArgs <- c('X', 'Y', 'formula', 'data')
  # mc[expectedArgs] <- lapply(mc[expectedArgs], eval.parent)
  
  ##============================= if data
  if (!is.null(try(mc$data)
  )) {
    ## ensure it's MAE class
    if (class(try(mc$data)
    )  !=  "MultiAssayExperiment") {
      .inv_data(mc$data)
    }
    
    ##--------------- if data & formula≠NULL
    ##--- i) if (data,formula) given change it to X and Y matrices
    if (class(try(mc$formula)
    )  !=  "NULL") {
      .sformula_checker(mc = mc)
      mc[c("X", "Y")] <- as.character(as.list(mc$formula)[3:2])
      mc <- .get_xy(mc, mcc)
    }
    ##--------------- if data & formula=NULL
    else {
      ## check X and Y exist
      if(any(sapply(mc[c("X", "Y")], function(xy) {class(try(xy))=="NULL"})))
        .inv_assay()
      ## in case they're stored in variables
      mc[c("X", "Y")] <- lapply( mc[c("X", "Y")], eval.parent)
      ## ensure it is a single character
      if(any(sapply( mc[c("X", "Y")], length)!=1))
        .stop(.subclass = "inv_XY", message = "'X' and 'Y' must be assay names from 'data'")
      mc <- .get_xy(mc, mcc)
    }
    ##--- if data, X and Y , expect X and Y to be assays and change them to matrices
    # else if(class(try(mc$formula))!="NULL"){
    #   ## if formula not a fomrula class, expect it to be NULL and X and Y to be assay/colData
    #   if(class(try(mc$formula))!="NULL")
    #     .stop(.subclass = "inv_formula", message = "'formula' must be a formula object of form Y~X")
    #
    # }
    
  }
  ##============================= if data=NULL and formula≠NULL
  else if (class(try(mc$formula))!="NULL"){
    mc$formula <- as.formula(mc$formula)
    .sformula_checker(mc=mc)
    mc[c('Y','X')] <- as.list(mc$formula)[2:3]
  }
  mc[c("X", "Y")] <- lapply( mc[c("X", "Y")], eval.parent)
  mc$data <- mc$formula <- NULL
  return(mc)
}