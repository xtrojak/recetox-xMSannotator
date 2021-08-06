num_nodes <<- 16
load("testdata/adduct_weights.rda")
adduct_weights <<- adduct_weights

patrick::with_parameters_test_that("Advanced annotation works:", {
  #skip("Currently excluded!")
  
  testname <<- test_identifier
  max.rt.diff <<- max_rt_diff

  wd <- getwd()

  peaks_filename <- paste0(testname, ".rds")
  peaks_filepath <- file.path(wd, "testdata", peaks_filename)

  outloc <- file.path(tempdir(), testname)

  peaks <- unique(readRDS(peaks_filepath))

  annotation <- multilevelannotation(
    peaks,
    num_nodes = num_nodes,
    outloc = outloc,
    db_name = database,
    queryadductlist = queryadductlist,
    adduct_weights = adduct_weights,
    max.rt.diff = max.rt.diff,
    allsteps = TRUE,
    cormethod = correlation_method,
    mass_defect_mode = mass_defect_mode,
    mode = mode
  )

  for (i in seq.int(from = 1, to = 5, by = 1)) {
    filename <- paste0("Stage", i, ".csv")
    actual <- read.csv(file.path(outloc, filename))
    expected <- read.csv(file.path(wd, "testdata", "advanced", testname, filename))

    actual <- dplyr::arrange(
      actual, mz, time,
    )
    expected <- dplyr::arrange(
      expected, mz, time
    )

    expect_equal(actual, expected)
  }

  setwd(wd)
},
cases(
    # qc_solvent = list(
    #   test_identifier = "qc_solvent",
    #   max_rt_diff = 0.5,
    #   queryadductlist = c("M+H", "M+2H", "M+H+NH4", "M+ACN+2H"),
    #   database = "HMDB",
    #   correlation_method = "pearson",
    #   mass_defect_mode = "pos",
    #   mode = "pos"
    # )
    # qc_matrix = list(
    #   test_identifier = "qc_matrix",
    #   max_rt_diff = 0.5,
    #   queryadductlist = c("M+H", "M+2H", "M+ACN+Na", "M+2ACN+H", "2M+H"),
    #   database = "KEGG",
    #   correlation_method = "spearman",
    #   mass_defect_mode = "both",
    #   mode = "pos"
    # ),
    # batch1_neg_hmdb = list(
    #   test_identifier = "batch1_neg",
    #   max_rt_diff = 0.5,
    #   queryadductlist = c("M-H", "M-2H"),
    #   database = "HMDB",
    #   correlation_method = "pearson",
    #   mass_defect_mode = "both",
    #   mode = "neg"
    # ),
    sample_data_custom = list(
      test_identifier = "sample_data_custom",
      max_rt_diff = 2,
      queryadductlist = c(
        "M+H", "M+2H", "M+H+NH4", "M+ACN+2H",
        "M+2ACN+2H", "M+NH4", "M+Na", "M+ACN+H",
        "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na",
        "2M+ACN+H", "M+2Na-H", "M+H-H2O", "M+H-2H2O"
      ),
      database = "HMDB",
      correlation_method = "pearson",
      mass_defect_mode = "pos",
      mode = "pos"
    )
  )
)