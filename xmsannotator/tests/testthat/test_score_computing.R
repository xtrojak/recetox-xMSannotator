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

patrick::with_parameters_test_that("abundance_ratio_computing", {
  actual <- compute_abundance_ratio(test_formula)
  expect_equal(actual, expected, tolerance = expected * 0.001)
},
  cases(hydrogen = list(
        test_formula = "H",
        expected = 0.00015
        ),
        chlorophorm = list(
        test_formula = "CHCl3",
        expected = 0.4126683
        ),
        chlorophyll = list(
        test_formula = "C55H72O5N4Mg",
        expected = 0.3171082
        )
  )
)