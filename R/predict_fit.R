#' Prediction function
#' 
#' @param fit An object of \code{stanfit}
#' @param data A data set
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
  D = matrix(NA, ncol = fitDATA$n_met, nrow = len_MCMC)
  for(m in 1:fitDATA$n_met){
    D[,m] =  fitMCMC$kem[,m] - (E + M) ;
  }
  # ACCUMULATION PHASE (0 =< t =< tacc)
  C0 <- fitDATA$C0
   
  conc <- matrix(NA, ncol = lentp, nrow = len_MCMC)
  if(n_met > 0){
    met <-  rep(NA, len_MCMC*lentp*n_met)
    dim(met) <- c(len_MCMC,lentp,n_met)
  } else{
    met <- NA
  }
  
   
  km <- fitMCMC$km
  kem <- fitMCMC$kem
  for(t in 1:rankacc){
    # Parent compound
    conc[,t] = (C0 - R[,t]) * exp(-(E + M) * time[t]) + R[t]
    # Metabolites
    if(n_met > 0){
      for(m in 1:n_met){
        met[,t,m] = km[,m] * (
          (C0-R[t])/ D[,m] * (exp(-(E + M)*time[t])-exp(-kem[,m] * time[t])) +
          R[t] / kem[,m] * (1 - exp(-(kem[,m] * time[t]))) 
        )
      }
    }
  }
  # DEPURATION PHASE (t > tacc)
  for(t in (rankacc+1):lentp){
    # Parent compound
    conc[,t] = (C0 - R[t] * (1 - exp((E + M)*tacc))) * exp(-(E + M) * time[t]) ;
    # Metabolites
    if(n_met > 0){
      for(m in 1:n_met){
        met[,t,m] = km[m] * (
          (C0-R[t]) / D[m] * (exp(-(E + M) * time[t]) - exp(-kem[m] * time[t])) + 
          R[t] / kem[m] * (exp(-kem[m] * (time[t]-tacc)) - exp(-kem[m] * time[t])) +
          R[t] / D[m] * (exp(-(E+M)*(time[t]-tacc)) - exp(-kem[m] * (time[t] - tacc)))
        )
      }
    }
  }
  # GROWTH
  if(n_out == 2){
    growth <- matrix(NA,  nrow = len_MCMC, ncol = lentp)
    G0 <- fitMCMC$G0
    gmax <- fitMCMC$gmax
    keg <- fitMCMC$ke[,2]
    for(t in 1:lentp){
      growth[,t] = (G0 - gmax) * exp(-keg * time[t]) + gmax
    }
  } else{
    growth <- NA
  }
  
  return(list(
    conc = conc,
    met = met,
    growth = growth
  ))
}
