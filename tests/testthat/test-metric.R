test_that("bioacc_metric", {
  
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
  
  BFC_MGS = bioacc_metric(fit_MGS, Male_Gammarus_Single)
  BFC_MGSG = bioacc_metric(fit_MGSG, Male_Gammarus_seanine_growth)
  BFC_CC = bioacc_metric(fit_CC, Chiro_Creuzot)
  
  # Check class
  expect_true(all(class(BFC_MGS) == c("bioaccMetric", "data.frame")))
  expect_true(all(class(BFC_MGSG) == c("bioaccMetric", "data.frame")))
  expect_true(all(class(BFC_CC) == c("bioaccMetric", "data.frame")))
  
  # Check colnames
  expect_true(all(colnames(BFC_MGS) == c("w")))
  expect_true(all(colnames(BFC_MGSG) == c("w")))
  expect_true(all(colnames(BFC_CC) == c("w","s","pw")))

  # Check class
  expect_true(all(class(plot(BFC_MGS)) == c("gg", "ggplot")))
  expect_true(all(class(plot(BFC_MGSG)) == c("gg", "ggplot")))
  expect_true(all(class(plot(BFC_CC)) == c("gg", "ggplot")))
  
})