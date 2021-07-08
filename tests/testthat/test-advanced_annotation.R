test_that("advanced annotation works", {
  tmpdir <- tempdir()
  load("testdata/sample_peaks.rda")
  load("testdata/expected_stage_1.rda")

  queryadductlist <- c(
    "M+H", "M+2H", "M+H+NH4", "M+ACN+2H",
    "M+2ACN+2H", "M+NH4", "M+Na", "M+ACN+H",
    "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na",
    "2M+ACN+H", "M+2Na-H", "M+H-H2O", "M+H-2H2O"
  )

  annotation <- multilevelannotation(
    peaks,
    num_nodes = 8,
    outloc = tmpdir,
    allsteps = FALSE,
    db_name = "HMDB",
    queryadductlist = queryadductlist,
    num_sets = 300,
  )

  actual_stage_1 <- read.csv("Stage1.csv")

  actual_stage_1 <- dplyr::arrange(actual_stage_1, Module_RTclust, mz, time)
  expected_stage_1 <- dplyr::arrange(expected_stage_1, Module_RTclust, mz, time)
  expect_equal(actual_stage_1, expected_stage_1)
})