test_that("test complex is running", {
  
  data("Male_Gammarus_Single")
  
  object = Male_Gammarus_Single
  time_accumulation = 4
  modelData_MGS = modelData(Male_Gammarus_Single, time_accumulation = 4)
  fit <- fitTK(modelData_MGS, iter = 100)
  print(fit)
  
  library(rstan)
  fitMCMC = rstan::extract(fit)
  quantile(fitMCMC$ku, probs = c(0.025, 0.5, 0.975))
  quantile(fitMCMC$ke, probs = c(0.025, 0.5, 0.975))
  quantile(fitMCMC$sigmaCpred, probs = c(0.025, 0.5, 0.975))
  
})