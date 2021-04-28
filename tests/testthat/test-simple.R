test_that("simple example", {
  
  library(rbioacc)
  data("Male_Gammarus_Single")
  test_object = Male_Gammarus_Single
  
  test_that("TEST .check_modelData_object", {
    expect_null(rbioacc:::.check_modelData_object(test_object))
  })
  test_that("TEST .modelDataSingle", {
    test_ls_object <- base::split(test_object, test_object$replicate)
    test_sgl_object = test_ls_object[[1]]
    test_ls_class = rbioacc:::.modelDataSingle(test_sgl_object, time_accumulation = 4)
    expect_equal(class(test_ls_class), "list")
  })
  
  modelData_MGS = modelData(test_object, time_accumulation = 4)
  
  stanData_MGS = modelData_MGS[[1]]
  
  stanData_MGS$Gobs = rep(0,stanData_MGS$lentp)
  
  fit <- complexTK_stan(stanData_MGS, iter = 25000)
  print(fit)
  
  library(rstan)
  fitMCMC = extract(fit)
  quantile(fitMCMC$ku, probs = c(0.025, 0.5, 0.975))
  quantile(fitMCMC$ke, probs = c(0.025, 0.5, 0.975))
  quantile(fitMCMC$sigmaCpred, probs = c(0.025, 0.5, 0.975))

})
