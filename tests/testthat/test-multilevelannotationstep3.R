patrick::with_parameters_test_that(
  "multilevelannotation step 3 works",
  {
    testthat_wd <- getwd()
    load("testdata/adduct_weights.rda")
    
    outloc <- file.path(tempdir(), subfolder)
    testdata_dir <- file.path(getwd(), "testdata", subfolder)
    # load expected results
    expected <- read.csv(file.path(testdata_dir, "Stage3.csv"))
    expected$MatchCategory <- as.character(expected$MatchCategory)

    # load needed objects
    chemscoremat <- readRDS(file.path(testdata_dir,"chemscoremat.Rds"))
    
    dir.create(outloc, recursive = TRUE)
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
      outloc=outloc,
       chemscoremat=chemscoremat, 
        adduct_weights=adduct_weights,
        pathwaycheckmode="pm"
    )

    actual <- read.csv(file.path(outloc, "Stage3.csv"))
    actual$MatchCategory <- as.character(actual$MatchCategory)
    #row.names(actual) <- 1:nrow(actual)

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
    
    # Annihilate
    setwd(testthat_wd)
    
    # compare with expected result
    expect_equal(actual, expected)
  },
  patrick::cases(
    qc_solvent = list(subfolder = "qc_solvent"),
    qc_matrix = list(subfolder = "qc_matrix"),
    batch1_neg = list(subfolder = "batch1_neg"),
    sourceforge = list(subfolder = "sourceforge")
  )
)