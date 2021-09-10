patrick::with_parameters_test_that(
  "multilevelannotation step 3 works",
  {
    testthat_wd <- getwd()
    testdata_dir <- file.path(testthat_wd, "testdata", subfolder)
    load(file.path(testdata_dir, "tempobjects.Rda"))

    outloc <- file.path(
      tempdir(),
      "multilevelannotationstep3",
      subfolder
    )
    dir.create(outloc, recursive = TRUE)

    # load expected results
    expected <- read.csv(file.path(testdata_dir, "Stage3.csv"))

    # load needed objects
    chemscoremat <- readRDS(file.path(testdata_dir, "chemscoremat.Rds"))

    file.copy(
      file.path(testdata_dir, "step1_results.Rda"),
      file.path(outloc, "step1_results.Rda")
    )
    file.copy(
      file.path(testdata_dir, "chemCompMZ.Rda"),
      file.path(outloc, "chemCompMZ.Rda")
    )

    # compute annotation step 3
    actual <- multilevelannotationstep3(
      outloc = outloc,
      chemscoremat = chemscoremat,
      adduct_weights = adduct_weights,
      pathwaycheckmode = pathwaycheckmode,
      num_sets=num_sets,
      boostIDs=boostIDs,
    )

    actual <- read.csv(file.path(outloc, "Stage3.csv"))

    setwd(testthat_wd)

    actual <- dplyr::arrange_all(actual)
    expected <- dplyr::arrange_all(expected)

    comparison <- dataCompareR::rCompare(
      actual,
      expected,
      keys = names(actual)
    )

    dataCompareR::saveReport(
      comparison,
      reportName = subfolder,
      reportLocation = outloc,
      showInViewer = FALSE,
      mismatchCount = 1000
    )

    expect_equal(actual, expected)
  },
  patrick::cases(
    qc_solvent = list(subfolder = "qc_solvent"),
    qc_matrix = list(subfolder = "qc_matrix"),
    batch1_neg = list(subfolder = "batch1_neg"),
    sourceforge = list(subfolder = "sourceforge")
  )
)