test_that("advanced annotation works", {
  tmpdir <- tempdir()
  load("testdata/sample_peaks.rda")
  load("testdata/expected_stage_1.rda")
  #load("testdata/expected_stage_2.rda")

  annotation <- multilevelannotation(
    peaks,
    num_nodes = 8,
    outloc = tmpdir,
    allsteps = FALSE,
    db_name = "HMDB"
  )

  actual_stage_1 <- read.csv("Stage1.csv")
  #actual_stage_2 <- read.csv("Stage2.csv")

  expect_equal(actual_stage_1, expected_stage_1)
  #expect_equal(actual_stage_2, expected_stage_2)
})