patrick::with_parameters_test_that("Advanced annotation works:", {
  #skip("Currently excluded!")

  num_nodes <- 16
  load("testdata/adduct_weights.rda")
  wd <- getwd()

  peaks_filename <- paste0(test_identifier, ".rds")
  peaks_filepath <- file.path(wd, "testdata", peaks_filename)

  outloc <- file.path(
    tempdir(),
    "advanced",
    test_identifier
  )
  dir.create(outloc, recursive = TRUE)

  peaks <- readRDS(peaks_filepath)
  peaks <- unique(peaks)

  annotation <- multilevelannotation(
    peaks,
    num_nodes = num_nodes,
    outloc = outloc,
    db_name = database,
    queryadductlist = queryadductlist,
    adduct_weights = adduct_weights,
    max.rt.diff = max_rt_diff,
    allsteps = TRUE,
    cormethod = correlation_method,
    mass_defect_mode = mass_defect_mode,
    mode = mode
  )

  setwd(wd)

  key_columns <- c('mz', 'time', 'Name', 'Adduct', 'Formula', 'chemical_ID', 'score', 'Confidence')

  for (i in seq.int(from = 1, to = 5, by = 1)) {
    filename <- paste0("Stage", i, ".csv")
    actual <- read.csv(file.path(outloc, filename))
    expected <- read.csv(file.path(wd, "testdata", test_identifier, filename))

    keys <- key_columns[which(key_columns %in% colnames(actual))]

    actual <- dplyr::arrange_at(actual, keys)
    expected <- dplyr::arrange_at(expected, keys)

    comparison <- dataCompareR::rCompare(
      actual,
      expected,
      keys = keys,
      mismatches = 1000
    )

    dataCompareR::saveReport(
        comparison,
        reportName = paste0(test_identifier, "_Stage", i),
        reportLocation = outloc,
        showInViewer = FALSE,
        mismatchCount = 1000
    )

    expect_equal(actual, expected, label = filename)
  }
},
patrick::cases(
    qc_solvent = list(
      test_identifier = "qc_solvent",
      max_rt_diff = 0.5,
      queryadductlist = c("M+H", "M+2H", "M+H+NH4", "M+ACN+2H"),
      database = "HMDB",
      correlation_method = "pearson",
      mass_defect_mode = "pos",
      mode = "pos"
    ),
    qc_matrix = list(
      test_identifier = "qc_matrix",
      max_rt_diff = 0.5,
      queryadductlist = c("M+H", "M+2H", "M+ACN+Na", "M+2ACN+H", "2M+H"),
      database = "KEGG",
      correlation_method = "spearman",
      mass_defect_mode = "both",
      mode = "pos"
    ),
    batch1_neg_hmdb = list(
      test_identifier = "batch1_neg",
      max_rt_diff = 0.5,
      queryadductlist = c("M-H", "M-2H"),
      database = "HMDB",
      correlation_method = "pearson",
      mass_defect_mode = "both",
      mode = "neg"
    ),
    sourceforge = list(
      test_identifier = "sourceforge",
      max_rt_diff = 2,
      queryadductlist = c("M+H"),
      database = "HMDB",
      correlation_method = "pearson",
      mass_defect_mode = "pos",
      mode = "pos"
    )
  )
)