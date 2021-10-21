test_that("Utils data reading: peak table processing.", {
  peak_table <- readRDS("test-data/utils/peak_table.rds")
  peak_table <- select(peak_table, c("peak", "mz", "rt"))
  actual <- as_peak_table(peak_table, intensities = FALSE)

  expect_equal(actual, peak_table)
})

test_that("Utils data reading: bad peak table processing.", {
  peak_table <- readRDS("test-data/utils/peak_table.rds")

  expect_error(as_peak_table(peak_table, intensities = TRUE), "is\\.numeric\\(data\\[\\[.*\\]\\]\\).*not TRUE")
})