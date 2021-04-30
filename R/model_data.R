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
#' 
#' 
modelData.data.frame <- function(object, time_accumulation, ...){
  
  .check_modelData_object(object)
  
  obj_colname = base::colnames(object)
  rtrn_ls = list()
  ########################
  # 0. GENERAL TIME VECTOR
  rtrn_ls$tp <- sort(unique(object$time))
  rtrn_ls$lentp <- length(rtrn_ls$tp)
  
  # 1. HANDLE REPLICATE
  ls_object <- base::split(object, object$replicate)
  rtrn_ls$n_rep <- length(ls_object)
  
  ls_object <- lapply(ls_object, .readapt_time, time_reference = rtrn_ls$tp)
  
  # 2. Exposure routes
  col_exposure <- .index_col_exposure(object)
  rtrn_ls$n_exp <- length(col_exposure)
  
  col_metabolite <- .index_col_metabolite(object)
  rtrn_ls$n_met <- length(col_metabolite)

  ls_object <- lapply(
    ls_object, .modelDataSingle,
    list(col_exposure=col_exposure, col_metabolite=col_metabolite)
  )
  
  Cexp <- do.call("cbind", sapply(ls_object, `[`, 1))
  dim(Cexp) <- c(rtrn_ls$lentp,rtrn_ls$n_exp,rtrn_ls$n_rep)
  if(all(sapply(1:rtrn_ls$n_rep, function(i){ .is_equal_rmInf(Cexp[,,1], Cexp[,,i]) } ))){
    rtrn_ls$Cexp <- as.matrix(Cexp[,,1])
  } else{
    stop("Replicates must have same Exposure profile")
  }

  rtrn_ls$Cobs = do.call("cbind", sapply(ls_object, `[`, 2))
  
  Cmet = do.call("cbind", sapply(ls_object, `[`, 3))
  dim(Cmet) <- c(rtrn_ls$lentp,rtrn_ls$n_met,rtrn_ls$n_rep)
  rtrn_ls$Cmet <- Cmet
  
  rtrn_ls$Gobs = do.call("cbind", sapply(ls_object, `[`, 4))
  
  # 3. Accumulation time
  rtrn_ls$tacc = time_accumulation
  rtrn_ls$rankacc = match(time_accumulation, rtrn_ls$tp)

  # 4. Prediction
  nsimu <- 500 # number of time points to simulate the model
  vt <- base::seq(0.0000001, max(rtrn_ls$tp)-0.0000001, length.out = nsimu)
  vt <- base::sort(c(vt, rtrn_ls$tp))
  rtrn_ls$vtacc = vt[vt<=rtrn_ls$tacc]
  rtrn_ls$len_vtacc = length(rtrn_ls$vtacc)
  rtrn_ls$vtdep = vt[vt>rtrn_ls$tacc]
  rtrn_ls$len_vtdep = length(rtrn_ls$vtdep)
  
  # 4. Priors Conc
  rtrn_ls$unifMax = 5*max(object$conc, na.rm = TRUE)
  
  # 5. Priors Growth
  if("growth" %in% colnames(object)){
    rtrn_ls$n_out <- 2
    rtrn_ls$gmaxsup <- 3*max(na.omit(object$growth, na.rm = TRUE))
  } else{
    rtrn_ls$n_out <- 1
    rtrn_ls$gmaxsup <- 0
  }
  
  rtrn_ls$C0 = mean(object[object$time == 0, ]$conc, na.rm = TRUE)

  return(rtrn_ls)
}



.is_equal_rmInf <- function(x,y){ 
  ux = unique(x) ; uy = unique(y)
  ux = ux[ux != Inf] ; uy = uy[uy != Inf]
  return(all(ux == uy))
}


.index_col_exposure <- function(object){
  col_exp = base::match(c("expw", "exps", "expf", "exppw"), base::colnames(object))
  return(col_exp[!base::is.na(col_exp)])
}

.index_col_metabolite <- function(object){
  obj_colname <- base::colnames(object)
  col_conc <- obj_colname[sapply(obj_colname, function(x){ regexpr("conc", x) == TRUE})]
  col_conc <- match(col_conc[col_conc != "conc"], obj_colname)
  return(col_conc[!base::is.na(col_conc)])
}


.readapt_time <- function(object, time_reference){
  obj_colname = base::colnames(object)
  lentp=length(time_reference)
  
  time_loc <- !(time_reference %in% object$time)
  if(sum(time_loc) > 0){
    time_NA <- time_reference[time_loc]
    # 1. time column as first column
    object <- cbind(time = object$time, object[, obj_colname != "time"])
    # 2. rebuild object with addition NA lines
    ls_NA <- lapply(1:length(time_NA), function(i){c(time_NA[i],rep(NA,ncol(object)-1))})
    row_NA <- do.call("rbind", ls_NA)
    colnames(row_NA) <- colnames(object)
    
    object <- rbind(object, row_NA)
    object <- object[order(object[, "time"]), ,drop = TRUE]
  }
  
  # convert NA to Inf because Stan does not support NA!
  object[is.na(object)] <- Inf
  return(object)
}

.modelDataSingle <- function(object, ls_col, ...){
  
  col_exp = ls_col$col_exposure
  col_conc = ls_col$col_metabolite
  
  rtrn_ls = list()
  # 1. Exposure
  rtrn_ls$Cexp = as.matrix(object[, col_exp])
  # 2. Parent Conc
  rtrn_ls$Cobs = as.matrix(object$conc)
  # 3. Metabolites
  if(length(col_conc) > 0){
    rtrn_ls$Cmet = as.matrix( object[, col_conc])
  } else{
    rtrn_ls$Cmet = matrix(0, nrow = nrow(object), ncol = 0)
  }
  # 4. Growth
  if("growth" %in% colnames(object)){
    rtrn_ls$Gobs = object$growth
  } else{
    rtrn_ls$Gobs = rep(0, nrow(object))
  }
  return(rtrn_ls)
}


.check_modelData_object <- function(object){
  
  obj_colname <- base::colnames(object)
  
  # time x replicate
  if(!('time' %in% obj_colname)) stop("`time` not a column name.")
  
  if(!('replicate' %in% obj_colname)) stop("`replicate` not a column name.")
  
  if(!('conc' %in% obj_colname)) stop("`conc` not a column name.")
  
  if(!any(c("expw", "exps", "expf", "exppw") %in% obj_colname)) stop("no exposure routes provided.")
  
}


