test_that("multilevelannotaton step 4 works", {
  
  subfolder = "qc_solvent"
  
  # load data needed during step 4
  load("../../data/adduct_weights.rda")
  
  # load expected results
  expected <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage4.csv"))
  expected$MatchCategory <- as.character(expected$MatchCategory)
  
  # set correct working directory
  testthat_wd <- getwd()
  test_path <- file.path(getwd(), "testdata/multilevelannotationstep4", subfolder)
  setwd(test_path)
  
  # compute annotation step 4
  result <- multilevelannotationstep4(outloc=".", max.mz.diff=10, max.rt.diff=0.5, filter.by="M+H", adduct_weights=adduct_weights)
  row.names(result) <- 1:nrow(result)
  
  # Annihilate
  setwd(testthat_wd)
  
  # compare with expected result
  expect_equal(result, expected)
  cmp <- arsenal::comparedf(result, expected, by = names(result))
  print(summary(cmp))
})
