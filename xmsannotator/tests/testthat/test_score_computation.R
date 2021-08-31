test_that("isotope_matching", {
  peaks <- readRDS("test-data/score_computation/peaks.Rda")
  query <- readRDS("test-data/score_computation/query.Rda")
  expected <- readRDS("test-data/score_computation/isotopes.Rda")

  actual <- compute_isotopes(peaks = peaks,
                             rt_tolerance = 1,
                             mass_defect_tolerance = 0.1,
                             max_isp = 5,
                             query)

  expect_equal(actual, expected, ignore_attr = TRUE)
})

patrick::with_parameters_test_that("abundance_ratio_computing", {
  actual <- compute_abundance_ratio(test_formula)
  expect_equal(actual, expected, tolerance = 0.01)
},
  patrick::cases(hydrogen = list(
        test_formula = "H2",
        expected = 0
        ),
        chlorophorm = list(
        test_formula = "CHCl3",
        expected = 0.9599
        ),
        chlorophyll = list(
        test_formula = "C55H72O5N4Mg",
        expected = 0.6
        )
  )
)

test_that("assign_abundance_ratios", {
  isotopes <- readRDS("test-data/score_computation/molecules.rda")
  expected <- readRDS("test-data/score_computation/molecules_with_abundances.rda")
  actual <- assign_isotope_abundances(isotopes)

  expect_equal(actual, expected)
})