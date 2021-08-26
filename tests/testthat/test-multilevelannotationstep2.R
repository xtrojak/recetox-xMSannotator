test_that("multilevelannotationstep2:", {
    test_identifier = "base"
    #skip("Currently excluded!")
    # Arrange
    load("testdata/adduct_weights.rda")

    testthat_wd <- getwd()
    test_path <- file.path(getwd(), "testdata/multilevelannotationstep2", test_identifier)
    setwd(test_path)

    load(file = "chem_score1.Rda")
    expected <- curchemscoremat
    curchemscoremat <- NA

    load(file = "step1_results.Rda")
    load(file = "global_cor.Rda")

    # Act
    actual <- multilevelannotationstep2(
        ".",
        1,
        max.rt.diff = max.rt.diff,
        chemids_split = chemids_split,
        num_sets = num_sets,
        mchemdata = mchemdata,
        mass_defect_window = mass_defect_window,
        corthresh = corthresh,
        global_cor = global_cor,
        mzid = mzid,
        adduct_table = adduct_table,
        adduct_weights = adduct_weights,
        max_isp = max_isp,
        MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
        mass_defect_mode = mass_defect_mode,
        chemids = chemids,
        isop_res_md = isop_res_md,
        filter.by = filter.by
    )
    #load(file = "chem_score1.Rda")

    # Assert
    expect_equal(actual, expected)

    # Annihilate
    setwd(test_path)
    unlink("stage2", recursive = TRUE)
    setwd(testthat_wd)
})
# },
# patrick::cases(
#     base = list(test_identifier = "base")
#     #cor_mz_matrix = list(test_identifier = "cor_mz_matrix")
# ))
