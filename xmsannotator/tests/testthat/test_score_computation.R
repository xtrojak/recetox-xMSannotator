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

test_that("Chemscore computing: isotope matching by intensity.", {
  query <- readRDS("test-data/score_computation/query.Rda")
  isotopes <- readRDS("test-data/score_computation/isotopes.Rda")
  pattern <- readRDS("test-data/score_computation/adenine_pattern.rda")

  expected <- readRDS("test-data/score_computation/matched_isotopes.rda")
  actual <- match_isotopes_by_intensity(query = query,
                                        isotopes = isotopes,
                                        pattern = pattern,
                                        intensity_deviation_tolerance = 0.4)
  expect_equal(actual, expected)
})

patrick::with_parameters_test_that("abundance_ratio_computing", {
  actual <- compute_isotopic_pattern(formula, minAbund = 0.001)
  expect_equal(actual, expected)
},
  patrick::cases(carbon = list(
        formula = "C",
        expected = readRDS("test-data/score_computation/C_pattern.rda")
        ),
        chlorophorm = list(
        formula = "CHCl3",
        expected = readRDS("test-data/score_computation/CHCl3_pattern.rda")
        ),
        chlorophyll = list(
        formula = "C55H72O5N4Mg",
        expected = readRDS("test-data/score_computation/CHLA_pattern.rda")
        )
  )
)

test_that("assign_abundance_ratios", {
  isotopes <- readRDS("test-data/score_computation/all_isotopes.Rda")
  expected <- readRDS("test-data/score_computation/isotopes_with_abundances.rda")
  actual <- assign_isotope_abundances(isotopes)

  expect_equal(actual, expected)
})