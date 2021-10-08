init.chemscoremat <- function(chemscoremat) {
  if (is.na(chemscoremat)) {
    chemscoremat <- read.csv("Stage4.csv")
  }
  chemscoremat <- as.data.frame(chemscoremat)
  chemscoremat$mz <- as.numeric(chemscoremat$mz)
  chemscoremat
}

increase.multimatches.score <- function(multimatch_features, adduct_weights) {
  features_adducts <- which(multimatch_features$Adduct %in% adduct_weights[, 1])
  if (length(features_adducts) > 0) {
    multimatch_features$score[features_adducts] <- (multimatch_features$score[features_adducts]) * 100
  }
  multimatch_features
}

get.features <- function(mz, match) {
  feature_count <- table(mz)
  features <- switch(match,
    "multiple" = feature_count[which(feature_count > 1)],
    "unique" = feature_count[which(feature_count == 1)]
  )
  features <- names(features)
}

remove.tmp.files <- function(loc) {
  fname <- paste("Stage5_annotation_results", sep = "")
  unlink(fname)
  suppressWarnings(unlink("*.Rda"))

  for (file in c("/stage2/", "/stage3/", "/stage4/", "/stage5/")) {
    loc <- paste(loc, file, sep = "")
    unlink(loc, force = TRUE, recursive = TRUE)
  }

  for (file in c("step1_results.Rda", "plot.pdf", "Rplots.pdf", "Rplots.pdf")) {
    try(unlink(file), silent = TRUE)
  }

  return(1)
}

multilevelannotationstep5 <- function(outloc,
                                      adduct_weights = NA,
                                      db_name = "HMDB",
                                      chemscoremat = NA,
                                      num_nodes = 2) {
  setwd(outloc)
  curated_res <- init.chemscoremat(chemscoremat)

  if (is.na(adduct_weights)) {
    adduct_weights <- data.frame(
      Adduct = c("M+H", "M-H"),
      Weight = c(1, 1)
    )
  }

  curated_res <- curated_res[order(curated_res$Confidence,
    curated_res$chemical_ID,
    curated_res$score,
    curated_res$Adduct,
    decreasing = TRUE
  ), ]

  duplicated_features <- get.features(curated_res$mz, "multiple")

  cl <- makeSOCKcluster(num_nodes)
  bad_ind <- {}
  bad_ind <- foreach(mz = duplicated_features, .combine = rbind) %dopar% {
    multimatches_idx <- which(curated_res$mz %in% mz)
    multimatch_features <- increase.multimatches.score(curated_res[multimatches_idx, ], adduct_weights)
    max_score_idx <- which(multimatch_features$score == max(multimatch_features$score, na.rm = TRUE))
    multimatches_idx <- multimatches_idx[-max_score_idx]

    # For all features with the same 'mz' except those that have the highest scores
    # Change the score of all annotations of a given molecule to:
    # 'highest score among annotations of the given molecule * [num(annotations of the molecule) - 1] / num(annotations of the molecule)
    # The purpose of this part of code appears to be increasing scores of isotopes.
    for (feature in 1:nrow(multimatch_features)) {
      if (multimatch_features$score[feature] != max(multimatch_features$score)) {
        same_molecule_idx <- which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[feature])
        same_molecule_annotation <- curated_res[same_molecule_idx, ]
        num_annotations <- nrow(same_molecule_annotation)
        score <- (num_annotations - 1) * max(same_molecule_annotation$score) / num_annotations
        curated_res$score[same_molecule_idx] <- score
      }
    }
    return(multimatches_idx)
  }
  stopCluster(cl)

  if (length(bad_ind) > 0) {
    curated_res <- curated_res[-c(bad_ind), ]
  }

  unique_features <- get.features(curated_res$mz, "unique")
  curated_res$MatchCategory <- rep("Multiple", dim(curated_res)[1])
  curated_res$MatchCategory[which(curated_res$mz %in% unique_features)] <- "Unique"

  write.csv(curated_res, file = "Stage5.csv", row.names = FALSE)
  remove.tmp.files(outloc)
  return(curated_res)
}
