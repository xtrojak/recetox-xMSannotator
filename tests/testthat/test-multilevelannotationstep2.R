patrick::with_parameters_test_that("multilevelannotationstep2:", {
    #skip("Currently excluded!")
    # Arrange
    testthat_wd <- getwd()
    test_path <- file.path(getwd(), "testdata/multilevelannotationstep2", test_identifier)
    setwd(test_path)

    load(file = "chem_score1.Rda")
    expected <- curchemscoremat
    curchemscoremat <- NA

    # Act
    multilevelannotationstep2(".", 1)
    load(file = "chem_score1.Rda")

    # Assert
    expect_equal(curchemscoremat, expected)

    # Annihilate
    setwd(test_path)
    unlink("stage2", recursive = TRUE)
    setwd(testthat_wd)
},
cases(
    base = list(test_identifier = "base")
    #cor_mz_matrix = list(test_identifier = "cor_mz_matrix")
))
