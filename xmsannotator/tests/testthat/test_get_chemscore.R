patrick::with_parameters_test_that("Compute chemscore can be called isolated", {
    # Arrange
    testthat_wd <- getwd()
    test_path <- file.path(
      testthat_wd,
      "test-data",
      test_identifier
    )
    load(file = file.path(test_path, "tempobjects.Rda"))

    outloc <- file.path(tempdir(), "get_chemscore", test_identifier)

    expected <- readRDS(file.path(test_path, "chemscoremat.Rds"))
    if(!dir.exists(outloc)) {
        dir.create(outloc, recursive = TRUE)
    }

    annotation_file <- file.path("test-data", "get_chemscore",paste0(test_identifier, "_annotation.Rds"))
    isotopes <- readRDS(annotation_file)

    if(is.factor(isotopes$MonoisotopicMass)) {
      isotopes$MonoisotopicMass <- as.numeric(levels(isotopes$MonoisotopicMass))[isotopes$MonoisotopicMass]
    }

    load(file = file.path(test_path, "global_cor.Rda"))

    actual <- purrr::pmap_dfr(
        isotopes,
        ~get_chemscore(...,
          annotation = isotopes,
          mass_defect_window = 0.01,
          adduct_weights = adduct_weights,
          corthresh = 0.7,
          global_cor = global_cor,
          max_diff_rt = max_diff_rt,
          outlocorig = outloc
        )
    )

    actual <- unique(actual)

    actual$Formula <- gsub(actual$Formula, pattern = "_.*", replacement = "")
    keys <- c("mz", "time", "Name", "Adduct", "Formula", "chemical_ID", "cur_chem_score")

    actual$MonoisotopicMass[is.na(actual$MonoisotopicMass)] <- "-"
    actual$theoretical.mz[is.na(actual$theoretical.mz)] <- "-"

    actual$MonoisotopicMass <- as.character(actual$MonoisotopicMass)
    actual$theoretical.mz <- as.character(actual$theoretical.mz)
    expected$Module_RTclust <- as.character(expected$Module_RTclust)

    actual <- dplyr::arrange_at(
      actual, keys
    )
    expected <- dplyr::arrange_at(
      expected, keys
    )

    comparison <- dataCompareR::rCompare(actual, expected, keys = keys)
    dataCompareR::saveReport(
      comparison,
      reportName = test_identifier,
      reportLocation = outloc,
      showInViewer = FALSE,
      mismatchCount = 10000
    )

    write.csv(actual, file = file.path(outloc, "chemscoremat_actual.csv"))
    write.csv(expected, file = file.path(outloc, "chemscoremat_expected.csv"))

    setwd(testthat_wd)

    # Assert
    expect_equal(actual, expected)

    rm(actual)
    gc(reset = TRUE)
},
  patrick::cases(
    qc_solvent = list(test_identifier = "qc_solvent", max_diff_rt = 0.5),
    batch1_neg = list(test_identifier = "batch1_neg", max_diff_rt = 0.5),
    sourceforge = list(test_identifier = "sourceforge", max_diff_rt = 2),
    qc_matrix = list(test_identifier = "qc_matrix", max_diff_rt = 0.5)
  )
)
