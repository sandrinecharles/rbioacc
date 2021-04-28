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
  
})