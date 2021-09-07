patrick::with_parameters_test_that(
  "multilevelannotation step 3 works",
  {
    testthat_wd <- getwd()
    load("../../data/adduct_weights.rda")
    
    # load expected results
    expected <- read.csv(file.path(getwd(), "testdata/advanced", subfolder, "Stage3.csv"))
    expected$MatchCategory <- as.character(expected$MatchCategory)
    
    # set correct working directory
    test_directory <- file.path(getwd(), "testdata/multilevelannotationstep3", subfolder)
    setwd(test_directory)
    
    # load needed objects
    load("stage2.Rda")
    
    # compute annotation step 3
    result <- multilevelannotationstep3(outloc=".", chemscoremat=chemscoremat, 
                                        adduct_weights=adduct_weights, pathwaycheckmode="pm")
    row.names(result) <- 1:nrow(result)
    
    # Annihilate
    setwd(testthat_wd)
    
    # compare with expected result
    cmp <- arsenal::comparedf(result, expected, by=names(result))
    expect_equal(arsenal::n.diff.obs(cmp), 0)
  },
  patrick::cases(
    qc_solvent = list(subfolder <- "qc_solvent"),
    qc_matrix = list(subfolder <- "qc_matrix"),
    batch1_neg = list(subfolder <- "batch1_neg"),
    sourceforge = list(subfolder <- "sourceforge")
  )
)