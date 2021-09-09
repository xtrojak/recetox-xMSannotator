test_that("multilevelannotationstep5", {
    skip("Currently excluded!")
    load("testdata/multilevelannotationstep5/tempobjects.Rda")
    load("testdata/adduct_weights.rda")

    outloc <- tempdir()
    file.copy("testdata/multilevelannotationstep5/Stage4.csv", file.path(outloc, "Stage4.csv"))

    expected <- read.csv("testdata/multilevelannotationstep5/Stage5.csv")
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
})