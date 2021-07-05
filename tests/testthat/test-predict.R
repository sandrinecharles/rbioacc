test_that("predict", {
  
  # small test to run fast
  data("Male_Gammarus_Single")
  modelData_MGS <- modelData(Male_Gammarus_Single, time_accumulation = 4)
  fit_MGS <- fitTK(modelData_MGS, iter = 10, chains = 2)
  data("Male_Gammarus_seanine_growth")
  modelData_MGSG <- modelData(Male_Gammarus_seanine_growth, time_accumulation = 1.417)
  fit_MGSG <- fitTK(modelData_MGSG, iter = 10, chains=2)
  data("Chiro_Creuzot")
  Chiro_Creuzot <- Chiro_Creuzot[Chiro_Creuzot$replicate == 1,]
  modelData_CC <- modelData(Chiro_Creuzot, time_accumulation = 1.0)
  fit_CC <- fitTK(modelData_CC, iter = 10, chains=2)
  
  # SAME data for test
  # Data 4 prediction should respect the exposure routes
  data_4pred_MGS <- data.frame( time = 1:25, expw = 4e-5 )
  predict_MGS <- predict(fit_MGS, data_4pred_MGS)

  data_4pred_MGSG <- data.frame( time = 1:21, expw = 18 )
  predict_MGSG <- predict(fit_MGSG, data_4pred_MGSG)
  
  data_4pred_CC <- data.frame( time = 1:25, expw = 18, exps = 1200, exppw = 15 )
  predict_CC <- predict(fit_CC, data_4pred_CC)

  test_true(all(class(plot(predict_MGS)) == c("gg","ggplot")))
  test_true(all(class(plot(predict_MGSG)) == c("gg", "ggplot")))
  test_true(all(class(plot(predict_CC)) == c("gg", "ggplot")))
  
})

  