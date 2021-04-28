#' Create a list giving data and parameters to use in the model inference.
#' @param object An object of class \code{data.frame}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
#' @return A \code{list} with data and parameters require for model inference.
#' 
#' 
modelData <- function(object, ...){
  UseMethod("modelData")
}


#' @rdname modelData
#' 
#' @export
#' 
#' 
modelData.data.frame <- function(object, time_accumulation, ...){
  
  .check_modelData_object(object)
  
  ls_object <- base::split(object, object$replicate)
  
  rtrn_ls = lapply(ls_object, .modelDataSingle, time_accumulation)
  
  C0 = mean(object[object$time == 0, ]$conc)
  rtrn_ls = lapply(rtrn_ls, append, list(C0=C0))
  
  return(rtrn_ls)
}


.modelDataSingle <- function(object, time_accumulation, ...){
  
  #################################
  # WORK ONLY FOR SINGLE REPLICATE
  ##################################
  obj_colname = base::colnames(object)
  rtrn_ls = list()
  # 0. Time vector
  rtrn_ls$tp = object$time
  rtrn_ls$lentp = length(rtrn_ls$tp)
  
  # 1. Exposure routes
  col_exp = base::match(c("expw", "exps", "expf", "exppw"), obj_colname)
  col_exp = col_exp[!base::is.na(col_exp)]
  
  rtrn_ls$Cexp = as.matrix(object[, col_exp])
  rtrn_ls$n_exp = ncol(rtrn_ls$Cexp)

  # 2. Parent Conc
  rtrn_ls$Cobs = object[, "conc"][[1]]
  
  # 3. Metabolites
  col_conc <- obj_colname[sapply(obj_colname, function(x){ regexpr("conc", x) == TRUE})]
  col_conc[col_conc != "conc"]
  col_conc = match(col_conc[col_conc != "conc"], obj_colname)
  col_conc = col_conc[!base::is.na(col_conc)]
  rtrn_ls$n_met = length(col_conc)
  
  if(rtrn_ls$n_met > 0){
    rtrn_ls$Cmet = as.matrix( object[, col_conc])
  } else{
    rtrn_ls$Cmet = matrix(0, nrow = rtrn_ls$lentp, ncol = 0)
  }
  
  # 4. Growth
  if("growth" %in% colnames(object)){
    rtrn_ls$n_out = 2
    rtrn_ls$gmaxsup <- 3*max(na.omit(object$growth))
    rtrn_ls$Gobs = object$growth
  } else{
    rtrn_ls$n_out = 1
    rtrn_ls$gmaxsup = 0
    rtrn_ls$Gobs = rep(0,rtrn_ls$lentp)
  }
  
  # 5. Accumulation time
  rtrn_ls$tacc = time_accumulation
  rtrn_ls$rankacc = match(time_accumulation, object$time)

  rtrn_ls$unifMax = 5*max(object$conc)
  
  # 6. Prediction
  nsimu <- 500 # number of time points to simulate the model
  vt <- base::seq(0.0000001, max(rtrn_ls$tp)-0.0000001, length.out = nsimu)
  vt <- base::sort(c(vt, object$time))
  rtrn_ls$vtacc = vt[vt<=rtrn_ls$tacc]
  rtrn_ls$len_vtacc = length(rtrn_ls$vtacc)
  rtrn_ls$vtdep = vt[vt>rtrn_ls$tacc]
  rtrn_ls$len_vtdep = length(rtrn_ls$vtdep)
  
  return(rtrn_ls)
}


.check_modelData_object <- function(object){
  
  obj_colname <- base::colnames(object)
  
  # time x replicate
  if(!('time' %in% obj_colname)) return("`time` not a column name.")
  
  if(!('replicate' %in% obj_colname)) return("`replicate` not a column name.")
  
  if(!('conc' %in% obj_colname)) return("`conc` not a column name.")
  
  if(!any(c("expw", "exps", "expf", "exppw") %in% obj_colname)) return("no exposure routes provided.")
  
}


