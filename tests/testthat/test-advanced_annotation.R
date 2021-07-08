test_that("advanced annotation Stage1 works", {
  tmpdir <- tempdir()
  load("testdata/sample_peaks.rda")
  peaks <- unique(peaks)

  load("testdata/expected_stage_1.rda")

  queryadductlist <- c(
    "M+H", "M+2H", "M+H+NH4", "M+ACN+2H",
    "M+2ACN+2H", "M+NH4", "M+Na", "M+ACN+H",
    "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na",
    "2M+ACN+H", "M+2Na-H", "M+H-H2O", "M+H-2H2O"
  )
  comparison_columns <- c("mz", "time", "rep1", "rep2", "rep3")

  annotation <- multilevelannotation(
    peaks,
    num_nodes = 8,
    outloc = tmpdir,
    allsteps = FALSE,
    db_name = "HMDB",
    queryadductlist = queryadductlist
  )

  actual_stage_1 <- read.csv("Stage1.csv")

  actual_stage_1 <- subset(actual_stage_1, select = -Module_RTclust)
  expected_stage_1 <- subset(expected_stage_1, select = -Module_RTclust)

  actual_stage_1 <- dplyr::arrange(
    actual_stage_1, mz, time, rep1, rep2, rep3
  )
  expected_stage_1 <- dplyr::arrange(
    expected_stage_1, mz, time, rep1, rep2, rep3
  )

  expect_equal(actual_stage_1, expected_stage_1)
})