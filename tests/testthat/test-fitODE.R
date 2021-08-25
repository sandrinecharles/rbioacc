test_that("fit ode", {
  
  data("Exposure_Sialis_lutaria")
  data("Internal_Sialis_lutaria")
  df_exposure = Exposure_Sialis_lutaria
  df_exposure$value = df_exposure$Cwater
  df_internal = Internal_Sialis_lutaria
  df_internal$value = df_internal$Cinternal

  modeldata_SL <- modelData_ode(df_exposure, df_internal, time_accumulation = 2.170)
  fit_SL <- fitTKvar(modeldata_SL, iter = 100)

})


