test_that("modelData outputs", {

  data("Male_Gammarus_Single")
  modelData(Male_Gammarus_Single, time_accumulation = 4)
  
  data("Male_Gammarus_seanine_growth")
  modelData(Male_Gammarus_seanine_growth, time_accumulation = 1.417)
  
  data("Male_Gammarus_Merged")
  modelData(Male_Gammarus_Merged, time_accumulation = 4)
  
  data("Oncorhynchus_two")
  modelData(Oncorhynchus_two, time_accumulation = 49)
  
  data("Chironomus_benzoapyrene")
  modelData(Chironomus_benzoapyrene, time_accumulation = 3)
  
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
  
})

