test_that("Integration utils: annotation table reformating", {
  main_annotation_table <- readRDS("test-data/integration_utils/main_annotation_table.rds")
  supplementary_data <- readRDS("test-data/integration_utils/supplementary_data.rds")

  actual <- reformat_annotation_table(main_annotation_table, supplementary_data)
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