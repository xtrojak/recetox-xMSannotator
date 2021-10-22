patrick::with_parameters_test_that("load_compound_table_parquet:", {
  filename <- file.path("test-data", "utils", paste0(db_name, ".parquet"))
  actual <- load_compound_table_parquet(filename)

  expected <- arrow::read_parquet(
    filename,
    col_select = c("monoisotopic_mass", "molecular_formula", "compound", "name")
  )
  expect_equal(actual, expected)
},
  db_name = c("hmdb", "pubchem")
)

test_that("Utils data reading: peak table processing.", {
  peak_table <- readRDS("test-data/utils/peak_table.rds")
  actual <- as_peak_table(peak_table, intensities = FALSE)
  peak_table <- select(peak_table, c("peak", "mz", "rt"))

  expect_equal(actual, peak_table)
})

test_that("Utils data reading: bad peak table processing.", {
  peak_table <- readRDS("test-data/utils/peak_table.rds")

  expect_error(as_peak_table(peak_table, intensities = TRUE), "is\\.numeric\\(data\\[\\[.*\\]\\]\\).*not TRUE")
})
