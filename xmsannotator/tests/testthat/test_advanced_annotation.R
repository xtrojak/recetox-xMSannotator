test_that("basic advanced_annotation functionality", {
  peaks <- arrow::read_parquet("test-data/advanced_annotation/batch1_neg.parquet")
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
    intensity_deviation_tolerance = 0.2
  )
  expect_success(annotation)
})
