test_that("Integration utils: annotation table reformatting", {
  main_annotation_table <- readRDS("test-data/integration_utils/main_annotation_table.rds")
  supplementary_data <- readRDS("test-data/integration_utils/supplementary_data.rds")

  actual <- reformat_annotation_table(isotopes_table, supplementary_data)
  expected <- readRDS("test-data/integration_utils/master_annotation_table.rds")

  expect_equal(actual, expected)
})