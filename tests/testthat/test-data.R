test_that("data load no complain", {
  
  expect_silent(data("Male_Gammarus_Single"))
  expect_silent(data("Male_Gammarus_Merged"))
  expect_silent(data("Male_Gammarus_seanine_growth"))
  
  expect_silent(data("Oncorhynchus_two"))
  expect_silent(data("Chironomus_benzoapyrene"))
  expect_silent(data("Chironomus_Creuzot"))
  
})

test_that("data.frame class", {
  
  expect_true("data.frame" %in% class(Male_Gammarus_Single))
  expect_true("data.frame" %in% class(Male_Gammarus_seanine_growth))
  expect_true("data.frame" %in% class(Male_Gammarus_Merged))
  expect_true("data.frame" %in% class(Oncorhynchus_two))
  expect_true("data.frame" %in% class(Chironomus_benzoapyrene))
  expect_true("data.frame" %in% class(Chironomus_Creuzot))
  
})