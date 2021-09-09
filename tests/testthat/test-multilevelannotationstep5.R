patrick::with_parameters_test_that("multilevelannotationstep5", {
    load("testdata/adduct_weights.rda")

    testdata_dir <- file.path("testdata", subfolder)
    load(file.path(testdata_dir, "tempobjects.Rda"))
    
    testthat_wd <- getwd()

    outloc <- file.path(tempdir(), subfolder, "multilevelannotationstep5")
    dir.create(outloc, recursive = TRUE)
    file.copy(
        file.path(testdata_dir,"Stage4.csv"),
        file.path(outloc, "Stage4.csv")
    )
    expected <- read.csv(file.path(testdata_dir, "Stage5.csv"))

    
    multilevelannotationstep5(
        outloc = outloc,
        max.rt.diff = max_diff_rt,
        filter.by = filter.by,
        adduct_weights = adduct_weights,
        min_ions_perchem = min_ions_perchem,
        boostIDs = boostIDs,
        max_isp = max_isp,
        db_name = db_name,
        max.mz.diff = max.mz.diff
    )

    actual <- read.csv(file.path(outloc, "Stage5.csv"))

    expect_equal(actual, expected)
    setwd(testthat_wd)
}, patrick::cases(
    qc_solvent = list(subfolder = "qc_solvent"),
    qc_matrix = list(subfolder = "qc_matrix"),
    batch1_neg = list(subfolder = "batch1_neg"),
    sourceforge = list(subfolder = "sourceforge")
))