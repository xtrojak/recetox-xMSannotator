patrick::with_parameters_test_that(
  "multilevelannotation step 4 works",
  {
    # load data needed during step 4
    load("../../data/adduct_weights.rda")
    
    chemscoremat <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage3.csv"))

    # load expected results
    expected <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage4.csv"))
    expected$MatchCategory <- as.character(expected$MatchCategory)

    # compute annotation step 4
    result <- multilevelannotationstep4(
      outloc = ".", chemscoremat = chemscoremat, max.mz.diff = 10, max.rt.diff = max_rt_diff,
      filter.by = "M+H", adduct_weights = adduct_weights
    )
    row.names(result) <- 1:nrow(result)

    # compare with expected result
    expect_equal(result, expected)
  },
  patrick::cases(
    qc_solvent = list(subfolder = "qc_solvent", max_rt_diff = 0.5),
    qc_matrix = list(subfolder = "qc_matrix", max_rt_diff = 0.5),
    batch1_neg = list(subfolder = "batch1_neg", max_rt_diff = 0.5),
    sourceforge = list(subfolder = "sourceforge", max_rt_diff = 2)
  )
)