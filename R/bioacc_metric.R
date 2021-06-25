#' @export
#' 
#' @rdname bioacc_metric
#' 
bioacc_metric <- function(fit, ...){
  UseMethod("bioacc_metric")
}

#' Biaccumulation metrics
#' 
#' @rdname bioacc_metric
#' 
#' @export
#' 
bioacc_metric.fitTK <- function(fit){
  
  fitMCMC <- rstan::extract(fit[["stanfit"]])
  sum_ <- apply(fitMCMC$ke, 1, sum) + apply(fitMCMC$km, 1, sum)
  BCF_ku <- lapply(1:ncol(fitMCMC$ku), function(i) fitMCMC$ku[,i] / sum_)
  
  return(BCF_ku)
}

#' Rertiev ku names from object
#'  
#' @export
#' 
ku_names <- function(object){
  col_exposure <- .index_col_exposure(object)
  sub <- substring(names(object[,col_exposure]), first = 4)
  return(sub)
}



