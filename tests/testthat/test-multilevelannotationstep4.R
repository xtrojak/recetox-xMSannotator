patrick::with_parameters_test_that("multilevelannotaton step 4 works", {
  
  # load data needed during step 4
  load("../../data/adduct_table.rda")
  load("../../data/adduct_weights.rda")
  
  # load expected results
  expected <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage4.csv"))
  expected$MatchCategory <- as.character(expected$MatchCategory)
  
  # set correct working directory
  testthat_wd <- getwd()
  test_path <- file.path(getwd(), "testdata/multilevelannotationstep4", subfolder)
  setwd(test_path)
  
  # compute annotation step 4
  result <- multilevelannotationstep4(".")
  row.names(result) <- 1:nrow(result)
  
  # Annihilate
  setwd(testthat_wd)
  
  # compare with expected result
  expect_equal(result, expected)
},
cases(
  qc_solvent = list(subfolder = "qc_solvent"),
  sample_data_custom = list(subfolder = "sample_data_custom"),
  batch1_neg = list(subfolder = "batch1_neg"),
  qc_matrix = list(subfolder = "qc_matrix")
))
