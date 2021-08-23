test_that("simple_annotation functionality on sample data", {
  peaks <- arrow::read_parquet("test-data/simple_annotation_sample/peaks.parquet")
  DB <- arrow::read_parquet("test-data/simple_annotation_sample/compounds.parquet")
  expected <- arrow::read_parquet("test-data/simple_annotation_sample/expected.parquet")
  
  result <- simple_annotation(peaks, DB)
  
  expect_equal(result, expected)
})


test_that("simple_annotation functionality on real data", {
  peaks <- arrow::read_parquet("test-data/simple_annotation_real/peaks.parquet")
  DB <- arrow::read_parquet("test-data/simple_annotation_real/compounds.parquet")
  adducts <- arrow::read_parquet("test-data/simple_annotation_real/adducts.parquet")
  expected <- arrow::read_parquet("test-data/simple_annotation_real/expected.parquet")

  result <- simple_annotation(peaks, DB, adduct_table=adducts, mass_tolerance=10e-6)
  
  names(result) <- c("Peak", "Adduct", "compound", "expected_mass", "mz", "time", "Name", "MonoisotopicMass", "Formula", "MatchCategory")
  comparison_columns = c("Adduct","Name", "mz", "time", "Formula", "MonoisotopicMass")
  result <- unique(subset(result, select = comparison_columns))
  
  cmp <- arsenal::comparedf(result, expected, by = c("Name", "mz", "time", "Formula", "Adduct", "MonoisotopicMass"))
  # check if number of differences is smaller than 4 - this is caused by typos/errors in reference data (produced by original version) 
  # and/or correct computation of check_element function (some formulas not containing C are misclassified, e.g. when containing Cl)
  expect_lte(arsenal::n.diffs(cmp), 4)
})
