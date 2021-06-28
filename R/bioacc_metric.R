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
bioacc_metric.fitTK <- function(fit, object){
  
  fitMCMC <- rstan::extract(fit[["stanfit"]])
  
  if(!is.null(fitMCMC$km)){
    sum_ <- apply(fitMCMC$ke, 1, sum) + apply(fitMCMC$km, 1, sum)
  } else{
    sum_ <- apply(fitMCMC$ke, 1, sum)
  }
  
  BCF_ku <- lapply(1:ncol(fitMCMC$ku), function(i) fitMCMC$ku[,i] / sum_)
  names(BCF_ku) <- exposure_names(object)
  
  df <- data.frame(do.call("cbind", BCF_ku))
  class(df) <- append("bioaccMetric", class(df))
  
  return(df)
}

#' Plot function for object of class \code{bioaccMetric}
#'  
#' @export
#' 
plot.bioaccMetric <- function(df){

  df_plt <- .fonte(df, "exposure","BCF_ku")
  
  df_quant <- .fonte(.df_quant95(df), "Quantile", "BCF_ku")
  df_quant$exposure <- rep(colnames(df), 3)
  
  plt <- ggplot() + 
    theme_minimal() +
    labs(x = "BCF", y = "Density") +
    geom_density(data = df_plt,
                 aes(x = BCF_ku, fill = exposure), fill = "grey", color = NA) +
    geom_vline(data = df_quant,
               aes(xintercept = BCF_ku, group = Quantile), linetype = "dashed") +
    facet_wrap(~exposure, scale = "free")
    
  return(plt)
}

#' Retrieve exposure routes names from object
#'  
#' @export
#' 
exposure_names <- function(object){
  col_exposure <- .index_col_exposure(object)
  sub <- substring(names(object[,col_exposure]), first = 4)
  return(sub)
}



# ------------------- INTERNAL

.fonte <- function(df, names_to, values_to){
  dfout <- data.frame(
    names_to = rep(colnames(df), each = nrow(df)),
    values_to = do.call("c", lapply(1:ncol(df), function(i) df[,i]))
  )
  colnames(dfout) <- c(names_to, values_to)
  return(dfout)
}
