test_that("isotope_matching", {
  peaks <- readRDS("testdata/peaks.Rda")
  query <- readRDS("testdata/query.Rda")
  expected <- readRDS("testdata/isotopes.Rda")


  actual <- compute_isotopes(peaks = peaks,
                             rt_tolerance = 1,
                             mass_defect_tolerance = 0.1,
                             max_isp = 5,
                             query)

  expect_equal(actual, expected, ignore_attr = TRUE)
})