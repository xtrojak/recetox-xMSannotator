test_that("advanced annotation Stage1 works", {
  skip("Currently excluded!")
  tmpdir <- tempdir()
  load("testdata/sample_peaks.rda")
  peaks <- unique(peaks)

  load("testdata/expected_stage_1.rda")

  comparison_columns <- c("mz", "time", "rep1", "rep2", "rep3")

  annotation <- multilevelannotation(
    peaks,
    num_nodes = 16,
    outloc = tmpdir,
    db_name = "HMDB",
    allsteps = FALSE
  )

  actual <- read.csv("Stage1.csv")

  actual <- subset(actual, select = -Module_RTclust)
  expected <- subset(expected_stage_1, select = -Module_RTclust)

  actual <- dplyr::arrange(
    actual, mz, time, rep1, rep2, rep3
  )
  expected <- dplyr::arrange(
    expected, mz, time, rep1, rep2, rep3
  )

  expect_equal(actual, expected)
})

test_that("advanced annotation Stage2 works", {
  # skip("Currently excluded!")
  here <- getwd()
  tmpdir <- tempdir()
  peaks <- readRDS("testdata/qc_solvent.rda")
  peaks <- unique(peaks)
  peaks <- dplyr::rename(subset(peaks, select = -feature), time = rt)

  load("testdata/expected_stage_2.rda")

  queryadductlist <- c(
    "M+H", "M+2H", "M+H+NH4", "M+ACN+2H",
    "M+2ACN+2H", "M+NH4", "M+Na", "M+ACN+H",
    "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na",
    "2M+ACN+H", "M+2Na-H", "M+H-H2O", "M+H-2H2O"
  )

  data(adduct_table)
  load("testdata/adduct_weights.rda")
  adduct_weights <- as.data.frame(adduct_weights)
  ls(adduct_weights)
  data(adduct_table)
  data(adduct_weights)

  num_nodes <- 8

  annotation <- multilevelannotation(
    peaks,
    num_nodes = num_nodes,
    outloc = tmpdir,
    db_name = "HMDB",
    queryadductlist = queryadductlist,
    adduct_weights = adduct_weights,
    max.rt.diff = 0.5,
    allsteps = TRUE
  )

  for (i in seq(from = 1, to = 5, by = 1)) {
    filename <- paste0("Stage", i, ".csv")
    actual <- read.csv(file.path(tmpdir, filename))
    expected <- read.csv(file.path(here, "testdata/advanced", filename))

    actual <- dplyr::arrange(
      actual, mz, time,
    )
    expected <- dplyr::arrange(
      expected, mz, time
    )

    expect_equal(actual, expected)
  }
})