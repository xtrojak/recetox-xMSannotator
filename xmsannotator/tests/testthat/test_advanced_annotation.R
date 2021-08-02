test_that("basic advanced_annotation functionality", {
  skip("Currently excluded!")
  peaks <- load_peak_table_hdf('aplcms_small.h5', TRUE)
  df <- advanced_annotation(peaks, sample_compound_table)

  expect_success(df)
})
