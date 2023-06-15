test_that("simple_annotation functionality on sample data", {
  peaks <- arrow::read_parquet("test-data/simple_annotation_sample/peaks.parquet")
  DB <- arrow::read_parquet("test-data/simple_annotation_sample/compounds.parquet")
  expected <- arrow::read_parquet("test-data/simple_annotation_sample/expected.parquet")
  
  result <- simple_annotation(peaks, DB)
  
  expect_equal(result, expected)
})
