test_that("multilevelannotation step 3 works", {
  
  subfolder <- "qc_solvent"
  testthat_wd <- getwd()
  load("../../data/adduct_weights.rda")
  
  # load expected results
  expected <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage3.csv"))
  expected$MatchCategory <- as.character(expected$MatchCategory)
  
  # set correct working directory
  test_directory <- file.path(getwd(), "testdata/multilevelannotationstep3", subfolder)
  setwd(test_directory)
  
  # load needed objects
  load("chem_score1.Rda")
  
  # compute annotation step 3
  result <- multilevelannotationstep3(outloc=".", chemscoremat=chemscoremat, 
                                      adduct_weights=adduct_weights, pathwaycheckmode="pm")
  row.names(result) <- 1:nrow(result)
  
  # Annihilate
  setwd(testthat_wd)
  
  # compare with expected result
  cmp <- arsenal::comparedf(result, expected, by=names(result))
  expect_equal(arsenal::n.diff.obs(cmp), 0)
})