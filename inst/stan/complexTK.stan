data {
  // Time points
  int<lower=0> lentp ;
  vector[lentp] tp ;

  // Exposure profiles
  int<lower=0> n_exp ;
  matrix[lentp, n_exp] Cexp ;
  
  // Internal concentraion
  vector[lentp] Cobs ;
  
  // Metabolites
  int<lower=0> n_met ;
  matrix[lentp, n_met] Cmet ;
  
  // Growth
  int<lower=0> n_out ;
  
  vector[lentp] Gobs ;
  real<lower=0> gmaxsup ;
  
  //
  int<lower=0> rankacc ;
  real<lower=0> tacc ;
  real<lower=0> C0 ;
  
  real unifMax ;
  
  int<lower=0> len_vtacc ; // length(vtacc)
  int<lower=0> len_vtdep ; // length(vtdep)
  vector[len_vtacc] vtacc ;
  vector[len_vtdep] vtdep ;
}
parameters {
  vector[n_exp] log10ku ; // uptake
  vector[n_out] log10ke ;
  vector[n_met] log10km ;
  vector[n_met] log10kem ;
  
  real<lower=0> sigmaCpred ;
  vector<lower=0>[n_met] sigmaCmetpred ;
  real<lower=0> sigmaGpred ;
  
  real<lower=0> gmax ;
  real<lower=0> G0 ;
  
}
transformed parameters{
  // PRIORS
  vector<lower=0>[n_exp] ku ;
  vector<lower=0>[n_out] ke ;
  vector<lower=0>[n_met] km ;
  vector<lower=0>[n_met] kem ;
  vector[lentp] U ;
  real M ;
  real E ;
  vector[lentp] R ;
  vector[n_met] D ;
  // little hack merging Cpred and Gpred
  matrix[lentp, n_out] CGpred ; 
  matrix[lentp,n_met] Cmetpred ;

  for(i in 1:n_exp){
    ku[i] = 10 ^ log10ku[i] ;
  }
  for(i in 1:n_out){
    ke[i] = 10 ^ log10ke[i] ;
  }
  for(i in 1:n_met){
    km[i] = 10 ^ log10km[i] ;
    kem[i] = 10 ^ log10kem[i] ;
  }

  M = sum(km) ;
  E = sum(ke) ;
  for(t in 1:lentp){
    // real operator*(row_vector x, vector y)
    U[t] =  Cexp[t,1:n_exp] * ku ;
    R[t] =  U[t] / (E + M) ;
  }
  for(i in 1:n_met){
      D[i] =  kem[i] - (E + M) ;
  }
  // ACCUMULATION PHASE (0 =< t =< tacc)
  for(t in 1:rankacc){
    // Parent compound
    CGpred[t, 1] = (C0 - R[t]) * exp(-(E + M) * tp[t]) + R[t] ;
    // Metabolites
    for(i in 1:n_met){
       Cmetpred[t,i] = km[i] * (
         (C0-R[t])/ D[i] * (exp(-(E+ M)*tp[t])-exp(-kem[i] * tp[t])) + R[t] / kem[i] * (1 - exp(-(kem[i] * tp[t]))) 
       ) ;
    }
  }
  //DEPURATION PHASE (t > tacc)
  for(t in (rankacc+1):lentp){
    // Parent compound
    CGpred[t, 1] = (C0 - R[t] * (1 - exp(E + M))) * exp(-(E + M) * tp[t]) ;
    // Metabolites
    for(i in 1:n_met){
      Cmetpred[t,i] = km[i] * (
        (C0-R[t]) / D[i] * (exp(-(E + M) * tp[t]) - exp(-kem[i] * tp[t])) + 
        R[t] / kem[i] * (exp(-kem[i] * (tp[t]-tacc)) - exp(-kem[i] * tp[t])) +
        R[t] / D[i] * (exp(-(E+M)*(tp[t]-tacc)) - exp(-kem[i] * (tp[t] - tacc)))
      ) ;
    }
  }
  // GROWTH
  if(n_out == 2){
    for(t in 1:lentp){
      CGpred[t, 2] = (G0 - gmax) * exp(-ke[2] * tp[t]) + gmax ;
    }
  }
}
model {
  // PRIORS
  target += uniform_lpdf(log10ku | -5, 5) ;
  target += uniform_lpdf(log10ke | -5, 5) ;
  target += uniform_lpdf(log10km | -5, 5) ;
  target += uniform_lpdf(log10kem | -5, 5) ;
  target += uniform_lpdf(sigmaCpred | 0, unifMax) ;
  target += uniform_lpdf(sigmaCmetpred | 0, unifMax) ;
  if(n_out == 2){
     target +=  uniform_lpdf(sigmaGpred | 0, unifMax) ;
     target +=  uniform_lpdf(gmax | gmaxsup/6, gmaxsup) ;
     target +=  uniform_lpdf(G0 | 0, gmaxsup) ;
  }
  
  // ACCUMULATION PHASE (0 =< t =< tacc) #
  for(t in 1:rankacc){
    // Parent compound
    target += normal_lpdf(Cobs[t] | CGpred[t, 1], sigmaCpred) ;
    // Metabolites
    for(i in 1:n_met){
       target += normal_lpdf(Cmet[t,i] | Cmetpred[t], sigmaCmetpred[i]) ;
    }
  }
  //DEPURATION PHASE (t > tacc)
  for(t in (rankacc+1):lentp){
    // Parent compound
    target += normal_lpdf(Cobs[t] | CGpred[t, 1], sigmaCpred) ;
    // Metabolites
    for(i in 1:n_met){
       target += normal_lpdf(Cmet[t,i] | Cmetpred[t], sigmaCmetpred[i]) ;
    }
  }
  if(n_out == 2){
      for(t in 1:lentp){
        target += normal_lpdf(Gobs[t] | CGpred[t, 2], sigmaGpred) ;
      }
  }
}
/*
generated quantities {
  // 0 < t < tacc
  vector[len_vtacc] Cpredp ;
  vector[len_vtacc] Cobsp ;
  matrix[len_vtacc,n_met] Cmetpredp ;
  matrix[len_vtacc,n_met] Cmetp ;
  // t > tacc
  vector[len_vtdep] Cpredpdep ;
  vector[len_vtdep] Cobspdep ;
  matrix[len_vtdep,n_met] Cmetpredpdep ;
  matrix[len_vtdep,n_met] Cmetpdep ;
  // growth
  vector[(len_vtacc+len_vtdep)] Gpredp ;
  //
  for(t in 1:len_vtacc){
    // Parent compound
    Cpredp[t] = (C0 - R[t]) * exp(-(E + M) * vtacc[t]) + R[t] ;
    Cobsp[t] = normal_rng(Cpredp[t], sigmaCpred) ;
    // Metabolites
    for(i in 1:n_met){
       Cmetpredp[t,i] = km[i] * (
         (C0-R[t])/ D[i] * (exp(-(E+ M)*vtacc[t])-exp(-kem[i] * vtacc[t])) + R[t] / kem[i] * (1 - exp(-(kem[i] * vtacc[t]))) 
       ) ;
       Cmetp[t,i] = normal_rng(Cmetpredp[t,i], sigmaCmetpred[i]) ;
    }
  }
  //
  for(t in 1:len_vtdep){
    // Parent compound
    Cpredpdep[t] = (C0 - R[t] * (1 - exp(E + M))) * exp(-(E + M) * vtdep[t]) ;
    Cobspdep[t] = normal_rng(Cpredpdep[t], sigmaCpred) ;
    // Metabolites
    for(i in 1:n_met){
      Cmetpredpdep[t,i] = km[i] * (
        (C0-R[t]) / D[i] * (exp(-(E + M) * vtdep[t]) - exp(-kem[i] * vtdep[t])) + 
        R[t] / kem[i] * (exp(-kem[i] * (vtdep[t]-tacc)) - exp(-kem[i] * vtdep[t])) +
        R[t] / D[i] * (exp(-(E+M)*(vtdep[t]-tacc)) - exp(-kem[i] * (vtdep[t] - tacc)))
      ) ;
      Cmetpdep[t,i] = normal_rng(Cmetpredpdep[t,i], sigmaCmetpred[i]) ;
    }
  }
  if(n_out == 2){
    for(t in 1:len_vtacc){
      Gpredp[t] = (G0 - gmax) * exp(-ke[2] * vtacc[t]) + gmax ;
    }
    for(t in (len_vtacc+1):(len_vtacc+len_vtdep)){
      Gpredp[t] = (G0 - gmax) * exp(-ke[2] * vtdep[t]) + gmax ;
    }
  }
}
*/

