test_that("isotope_matching", {
  peaks <- readRDS("testdata/peaks.Rda")
  query <- readRDS("testdata/query.Rda")
  expected <- readRDS("testdata/isotopes.Rda")

  actual <- compute_isotopes(peaks = peaks, query = query,
                             rt_tolerance = 1, mass_defect_tolerance = 0.1,
                             abundance_ratio = 0.2, max_isp = 5)

  expect_equal(actual, expected, ignore_attr = TRUE)
})