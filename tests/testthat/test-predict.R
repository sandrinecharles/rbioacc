test_that("predict", {
  
  # small test to run fast
  data("Male_Gammarus_Single")
  modelData_MGS <- modelData(Male_Gammarus_Single, time_accumulation = 4)
  fit_MGS <- fitTK(modelData_MGS, iter = 10, chains = 4)
  data("Male_Gammarus_seanine_growth")
  modelData_MGSG <- modelData(Male_Gammarus_seanine_growth, time_accumulation = 1.417)
  fit_MGSG <- fitTK(modelData_MGSG, iter = 10, chains=4)
  data("Chiro_Creuzot")
  Chiro_Creuzot <- Chiro_Creuzot[Chiro_Creuzot$replicate == 1,]
  modelData_CC <- modelData(Chiro_Creuzot, time_accumulation = 1.0)
  fit_CC <- fitTK(modelData_CC, iter = 10, chains=4)
  
  # SAME data for test
  data_MGS <- fit_MGS$origin_data
  data_MGSG <- fit_MGSG$origin_data
  data_CC <- fit_CC$origin_data
  
  predict(fit_MGS, data_MGS)
  
})