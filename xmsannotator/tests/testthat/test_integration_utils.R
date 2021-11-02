test_that("Integration utils: annotation table reformating", {
  main_annotation_table <- readRDS("test-data/integration_utils/main_annotation_table.rds")

  actual <- reformat_annotation_table(main_annotation_table)
  expected <- readRDS("test-data/integration_utils/master_annotation_table.rds")

  actual <- dplyr::arrange_all(actual)
  expected <- dplyr::arrange_all(expected)

  comparison <- dataCompareR::rCompare(
    actual,
    expected,
    keys = names(actual)
  )

  expect_equal(actual, expected)
})

patrick::with_parameters_test_that("Global cor adapter works", {
  testthat_wd <- getwd()
  test_path <- file.path(testthat_wd, "test-data")

  load(file.path(test_path, test_identifier, "global_cor.Rda"))

  peaks_filename <- file.path(test_path, paste0(test_identifier, ".parquet"))
  peak_table <- load_peak_table_parquet(peaks_filename)
  peak_intensity_matrix <- get_peak_intensity_matrix(peak_table)
  peak_correlation_matrix <- compute_peak_correlations(peak_intensity_matrix)

  actual <- reformat_correlation_matrix(peak_table, peak_correlation_matrix, truncate = TRUE)

  expect_true(all.equal(actual, global_cor, check.names = TRUE, countEQ = TRUE))
},
  patrick::cases(
    qc_solvent = list(test_identifier = "qc_solvent"),
    batch1_neg = list(test_identifier = "batch1_neg"),
    sourceforge = list(test_identifier = "sourceforge"),
    qc_matrix = list(test_identifier = "qc_matrix")
  )
)
