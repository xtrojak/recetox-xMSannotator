load_expected <- function(test_path) {
  load(file = file.path(test_path, "step1_results.Rda"))
  load(file = file.path(test_path, "global_cor.Rda"))
  load(file = file.path(test_path, "tempobjects.Rda"))

  outloc <- file.path(tempdir(), "get_chemscore", basename(test_path))

  if(!dir.exists(outloc)) {
      dir.create(outloc, recursive = TRUE)
  }

  expected <- multilevelannotationstep2(
      outloc = outloc,
      max_diff_rt = max.rt.diff,
      mchemdata = mchemdata,
      mass_defect_window = mass_defect_window,
      corthresh = corthresh,
      global_cor = global_cor,
      adduct_weights = adduct_weights,
      max_isp = max_isp,
      MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
      mass_defect_mode = mass_defect_mode,
      chemids = chemids,
      isop_res_md = isop_res_md,
      filter.by = filter.by
    )

    return(expected)
}

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

    if(!dir.exists(outloc)) {
      dir.create(outloc, recursive = TRUE)
    }

    # expected <- readRDS(file.path(test_path, "chemscoremat.Rds"))
    expected <- unique(load_expected(test_path))

    setwd(testthat_wd)

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
          mass_defect_window = mass_defect_window,
          adduct_weights = adduct_weights,
          corthresh = corthresh,
          global_cor = global_cor,
          max_diff_rt = max_diff_rt,
          outlocorig = outloc
        )
    )
    actual <- unique(actual)

    #actual$Formula <- gsub(actual$Formula, pattern = "_.*", replacement = "")
    keys <- c("mz", "time", "Name", "Adduct", "Formula", "chemical_ID", "cur_chem_score", "MatchCategory")

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
