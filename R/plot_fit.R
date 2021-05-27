.df_quant95 <- function(x){
  mat_quants = apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975))
  df = data.frame(
    q50 = mat_quants[2,],
    qinf95 = mat_quants[1,],
    qsup95 = mat_quants[3,]
  )
  return(df)
}

.add_data = function(df_quant95,data,id){
  if(is.vector(data)){
    df_quant95$time = fit$stanTKdata$tp
    df_quant95$observation = data
    df_quant95$replicate = 1
    df <- df_quant95
  } else{
    ls <- lapply(1:ncol(data),
                 function(i){
                   df_quant95$time = fit$stanTKdata$tp
                   df_quant95$observation = data[,i]
                   df_quant95$replicate = i
                   return(df_quant95)
                 })
    df <- do.call("rbind", ls)
  }
  df <- df[df$observation != Inf,]
  df$variable <-  id
  return(df)
}

.df_for_plot <- function(fit){
  fitMCMC = rstan::extract(fit[["stanfit"]])
  # 
  ls_out <- list()
  ls_out$conc <- .add_data(
    .df_quant95(fitMCMC$CGobs_out[,,1]),
    fit$stanTKdata$CGobs[,1,],
    "conc"
  )
  if(dim(fitMCMC$CGobs_out)[3] == 2){
    ls_out$growth <- .add_data(
      .df_quant95(fitMCMC$CGobs_out[,,2]),
      fit$stanTKdata$CGobs[,2,],
      "growth"
    )
  }
  if("Cmet_out" %in% names(fitMCMC)){
    for(i in 1:fit$stanTKdata$n_met){
      name <- paste0("concm",i)
      ls_out[[name]] <- .add_data(
        .df_quant95(fitMCMC$Cmet_out[,,i]),
        fit$stanTKdata$Cmet[,i,],
        name
      )
    }
  }
  
  df <- do.call("rbind", ls_out)
  
  return(df)
}

#' Plotting method for \code{fitTK} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{fitTK}.  It plots the fit obtained for each
#' variable in the original dataset.
#' 
#' @param fit And object returned by fitTK
#' 
#' @export
#' 
#' @import ggplot2
#' 
plot_fitTK <- function(fit){
#plot.fitTK <- function(fit){
  
  df <- .df_for_plot(fit)
  
  plt <- ggplot(data = df) + 
    theme_classic() +
    labs(x = "Time", y = "Concentration") +
    geom_ribbon(
      aes(x = time, ymin = qinf95, ymax = qsup95), fill = "grey80") +
    geom_line(aes(x = time, y = q50), color = "orange") +
    geom_point(aes(x = time, y = observation )) + 
    facet_wrap(~variable, scales = "free")

  return(plt)
  
}

