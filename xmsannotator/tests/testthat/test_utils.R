patrick::with_parameters_test_that("load_compound_table_parquet:", {
  filename <- file.path("test-data", "utils", paste0(db_name, ".parquet"))
  actual <- load_compound_table_parquet(filename)

  expected <- arrow::read_parquet(filename)
  expect_equal(actual, expected)
},
  db_name = c("hmdb", "pubchem")
)