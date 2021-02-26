test_that("simple_annotation", {
  data("sample_peak_table")
  data("sample_adduct_table")
  data("sample_compound_table")

  df <- simple_annotation(
    sample_peak_table,
    sample_compound_table,
    sample_adduct_table
  )

  expect_s3_class(df, 'data.frame')
})
