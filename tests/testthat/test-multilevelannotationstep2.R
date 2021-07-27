test_that("multilevelannotationstep2", {
    # Arrange
    setwd("testdata/multilevelannotationstep2")
    load(file = "chem_score1.rda")
    expected <- curchemscoremat
    curchemscoremat <- NA

    # Act
    multilevelannotationstep2(".", 1)
    load(file = "chem_score1.Rda")

    # Assert
    expect_equal(curchemscoremat, expected)

    # Annihilate
    setwd("../../..")
    unlink("stage2", recursive = TRUE)
})
