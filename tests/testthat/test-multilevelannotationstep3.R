test_that("multilevelannotation step 3 works", {
  
  subfolder <- "qc_solvent"
  testthat_wd <- getwd()
  
  # load expected results
  expected <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage3.csv"))
  expected$MatchCategory <- as.character(expected$MatchCategory)
  
  # set correct working directory
  test_directory <- file.path(getwd(), "testdata/multilevelannotationstep3", subfolder)
  dir.create(test_directory, recursive = TRUE)
  
  list.of.files <- list.files(file.path(getwd(), "testdata/multilevelannotationstep2/base"))
  
  file.copy(file.path(getwd(), "testdata/multilevelannotationstep2/base", list.of.files), test_directory)
  setwd(test_directory)
  
  # load needed objects
  load("chem_score1.Rda")
  load("tempobjects.Rda")
  
  # compute annotation step 3
  result <- multilevelannotationstep3(outloc=".",adduct_weights=adduct_weights,
                                      boostIDs=boostIDs,pathwaycheckmode=pathwaycheckmode)
  row.names(result) <- 1:nrow(result)
  
  # Annihilate
  setwd(testthat_wd)
  unlink(file.path(getwd(), "testdata/multilevelannotationstep3"), recursive = TRUE)
  
  # compare with expected result
  expect_equal(result, expected)
})