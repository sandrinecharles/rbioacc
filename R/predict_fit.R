#' Prediction function
#' 
#' @param fit An object of \code{stanfit}
#' @param data A data set with one column \code{time} and 1 to 4 exposure 
#' columns with name in \code{expw}, \code{exps}, \code{expf} and \code{exppw}
#' 
#' @export
#' 
predict.fitTK <- function(fit, data){
  
  fitDATA <- fit[["stanTKdata"]]
  fitMCMC <- rstan::extract(fit[["stanfit"]])
  
  n_met <- fitDATA$n_met
  len_MCMC <- nrow(fitMCMC$ku)
  # Exposure match
  data_origin <- fitDATA$origin_data
  
  col_names_origin <- colnames(data_origin)[ .index_col_exposure(data_origin)]
  Cexp_origin <- as.data.frame(data_origin[, col_names_origin])
  colnames(Cexp_origin) <- col_names_origin
  
  col_names <- colnames(data)[ .index_col_exposure(data)]
  Cexp <- as.data.frame(data[, col_names])
  colnames(Cexp) <- col_names
  
  if(!(all(colnames(Cexp) %in% colnames(Cexp_origin)) && 
       all(colnames(Cexp_origin) %in% colnames(Cexp)))){
    break("Exposure routes differ between 'fit' and 'data'")
  }
  
  n_exp <- ncol(Cexp)
  n_out <- fitDATA$n_out

  lentp <- nrow(data)
  time <- data$time
  
  tacc <- fitDATA$tacc
  rankacc <- match(tacc, time)
  
  if(n_met == 0){
    M = 0
  } else{
    M = fitMCMC$M
  }
  E = fitMCMC$E
   
  U = matrix(NA, ncol = lentp, nrow = len_MCMC)
  R = matrix(NA, ncol = lentp, nrow = len_MCMC)
  for(t in 1:lentp){
    for(i in 1:len_MCMC){
      U[i,t] =  sum(Cexp[t,] * fitMCMC$ku[i,])
    }
    R[,t] =  U[,t] / (E + M) ;
  }
  if(n_met > 0){
    D = matrix(NA, ncol = fitDATA$n_met, nrow = len_MCMC)
    for(m in 1:n_met){
      D[,m] =  fitMCMC$kem[,m] - (E + M) ;
    }
  }
  # ACCUMULATION PHASE (0 =< t =< tacc)
  C0 <- fitDATA$C0
   
  CGobs_out <- rep(NA, len_MCMC*lentp*n_out)
  dim(CGobs_out) <- c(len_MCMC,lentp,n_out)
  
  if(n_met > 0){
    Cmet_out <-  rep(NA, len_MCMC*lentp*n_met)
    dim(Cmet_out) <- c(len_MCMC,lentp,n_met)
  } else{
    Cmet_out <- NA
  }
  
  km <- fitMCMC$km
  kem <- fitMCMC$kem
  for(t in 1:rankacc){
    # Parent compound
    CGobs_out[,t,1] = (C0 - R[,t]) * exp(-(E + M) * time[t]) + R[,t]
    # Metabolites
    if(n_met > 0){
      for(m in 1:n_met){
        Cmet_out[,t,m] = km[,m] * (
          (C0-R[t])/ D[,m] * (exp(-(E + M)*time[t])-exp(-kem[,m] * time[t])) +
          R[,t] / kem[,m] * (1 - exp(-(kem[,m] * time[t])))
        )
      }
    }
  }
  # DEPURATION PHASE (t > tacc)
  for(t in (rankacc+1):lentp){
    # Parent compound
    CGobs_out[,t,1] = (C0 - R[,t] * (1 - exp((E + M)*tacc))) * exp(-(E + M) * time[t]) ;
    # Metabolites
    if(n_met > 0){
      for(m in 1:n_met){
        Cmet_out[,t,m] = km[m] * (
          (C0-R[t]) / D[m] * (exp(-(E + M) * time[t]) - exp(-kem[m] * time[t])) + 
          R[,t] / kem[m] * (exp(-kem[m] * (time[t]-tacc)) - exp(-kem[m] * time[t])) +
          R[,t] / D[m] * (exp(-(E+M)*(time[t]-tacc)) - exp(-kem[m] * (time[t] - tacc)))
        )
      }
    }
  }
  # GROWTH
  if(n_out == 2){
    G0 <- fitMCMC$G0
    gmax <- fitMCMC$gmax
    keg <- fitMCMC$ke[,2]
    for(t in 1:lentp){
      CGobs_out[,t,2] = (G0 - gmax) * exp(-keg * time[t]) + gmax
    }
  }
  
  predict_out <- list(
    time = time,
    CGobs_out = CGobs_out,
    Cmet_out = Cmet_out
  )
  
  class(predict_out) <- append("predictTK", class(predict_out))
  
  return(predict_out)
}
