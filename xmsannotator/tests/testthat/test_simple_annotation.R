test_that("simple_annotation functionality on sample data", {
  peaks <- arrow::read_parquet("test-data/simple_annotation_sample/peaks.parquet")
  DB <- arrow::read_parquet("test-data/simple_annotation_sample/compounds.parquet")
  expected <- arrow::read_parquet("test-data/simple_annotation_sample/expected.parquet")
  
  result <- simple_annotation(peaks, DB)
  
  expect_equal(result, expected)
})


test_that("simple_annotation functionality on real data", {
  skip("Currently excluded due to long computational time.")
  peaks <- arrow::read_parquet("test-data/simple_annotation_real/peaks.parquet")
  DB <- arrow::read_parquet("test-data/simple_annotation_real/compounds.parquet")
  expected <- arrow::read_parquet("test-data/simple_annotation_real/expected.parquet")

  result <- simple_annotation(peaks, DB, mass_tolerance=10e-6)
  
  names(result) <- c("Peak", "Adduct", "compound", "expected_mass", "mz", "time", "MonoisotopicMass", "Formula", "MatchCategory")
  
  comparison_columns = c("Adduct","Name", "mz", "time", "Formula", "MonoisotopicMass")
  expected <- subset(expected, select = comparison_columns)
  actual <- subset(annotation, select = comparison_columns)

  expect_equal(actual, expected)
})
