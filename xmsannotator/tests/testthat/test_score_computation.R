patrick::with_parameters_test_that("Chemscore computing: isotope matching", {
  peaks <- readRDS("test-data/score_computation/peaks.rds")
  actual <- compute_isotopes(peaks = peaks,
                             intensity_deviation_tolerance = 0.25,
                             rt_tolerance = 10,
                             mass_defect_tolerance = 0.1,
                             query)
  expect_equal(actual, expected, ignore_attr = TRUE)
},
  patrick::cases(
    isoproturon = list(
      query = readRDS("test-data/score_computation/isoproturon.rds"),
      expected = readRDS("test-data/score_computation/isoproturon_isp.rds")
    ),
    desisopropylatrazine = list(
      query = readRDS("test-data/score_computation/desisopropylatrazine.rds"),
      expected = readRDS("test-data/score_computation/desisopropylatrazine_isp.rds")
    ),
    acetochlor = list(
      query = readRDS("test-data/score_computation/acetochlor.rds"),
      expected = readRDS("test-data/score_computation/acetochlor_isp.rds")
    )
  )
)

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

patrick::with_parameters_test_that("Chemscore computing: theoretical isotopic pattern", {
  actual <- compute_isotopic_pattern(formula, minAbund = 0.001)
  expect_equal(actual, expected)
},
  patrick::cases(
    carbon = list(
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
