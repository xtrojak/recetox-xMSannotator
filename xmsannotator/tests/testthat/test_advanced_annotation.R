test_that("basic advanced_annotation functionality", {
  skip_on_ci()
  
  outloc <- file.path(tempdir(), "advanced_annotation")
  dir.create(outloc, recursive = TRUE)

  expected <- read.csv("test-data/advanced_annotation/expected.csv")

  peaks <- arrow::read_parquet("test-data/batch1_neg.parquet")
  DB <- arrow::read_parquet("test-data/advanced_annotation/hmdb.parquet")
  adduct_names <- c("M-H", "2M-H", "M-2H")
  adduct_table <- sample_adduct_table %>% filter(adduct %in% adduct_names)

  annotation <- advanced_annotation(
    peaks,
    DB,
    adduct_table = adduct_table,
    peak_rt_width = 15,
    time_tolerance = 15,
    min_cluster_size = 10,
    intensity_deviation_tolerance = 0.2,
    outloc = outloc
  )

  actual <- read.csv(file.path(outloc, "Stage5.csv"))

  expect_equal(actual, expected)
})
