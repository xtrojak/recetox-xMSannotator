init.chemscoremat <- function(chemscoremat) {
  if (is.na(chemscoremat)) {
    chemscoremat <- read.csv("Stage4.csv")
  }
  chemscoremat <- as.data.frame(chemscoremat)
  chemscoremat$mz <- as.numeric(chemscoremat$mz)
  chemscoremat
}

reevaluate.multimatches.score <- function(multimatch_features, adduct_weights) {
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
  bad_ind <- foreach(mz_idx = 1:length(duplicated_features), .combine = rbind) %dopar% {
    mz <- duplicated_features[mz_idx]
    common_mz_idx <- which(curated_res$mz %in% mz)

    multimatch_features <- reevaluate.multimatches.score(curated_res[common_mz_idx, ], adduct_weights)
    good_idx <- which(multimatch_features$score == max(multimatch_features$score, na.rm = TRUE))

    # For all features with the same 'mz' except those that have the highest scores
    # Change the score of all annotations of a given molecule to:
    # 'highest score among annotations of the given molecule * [num(annotations of the molecule) - 1] / num(annotations of the molecule)
    for (com_indval in 1:length(common_mz_idx)) {
      if (!com_indval %in% good_idx) {
        alike_annotations <- curated_res[which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval]), ]
        scoreval <- ((dim(alike_annotations)[1]) - 1) * alike_annotations$score[1] / (dim(alike_annotations)[1])
        scorevec <- c(rep(scoreval, length(which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval]))))
        curated_res$score[which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval])] <- scorevec
      }
    }
    common_mz_idx <- common_mz_idx[-good_idx]

    return(common_mz_idx)
  }
  stopCluster(cl)

  if (length(bad_ind) > 0) {
    curated_res <- curated_res[-c(bad_ind), ]
  }

  unique_features <- get.features(curated_res$mz, "unique")

  curated_res$MatchCategory <- rep("Multiple", dim(curated_res)[1])
  curated_res$MatchCategory[which(curated_res$mz %in% unique_features)] <- "Unique"

  write.csv(curated_res, file = "Stage5.csv", row.names = FALSE)

  curated_res$mz <- sprintf("%.5f", as.numeric(curated_res$mz))





  chemIDs <- curated_res$chemical_ID
  htmllink <- switch(db_name,
    "HMDB" = paste("<a href=http://www.hmdb.ca/metabolites/", chemIDs, ">", chemIDs, "</a>", sep = ""),
    "KEGG" = paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", chemIDs, ">", chemIDs, "</a>", sep = ""),
    "LipidMaps" = paste("<a href=http://www.lipidmaps.org/data/LMSDRecord.php?LMID=", chemIDs, ">", chemIDs, "</a>", sep = ""),
    "T3DB" = paste("<a href=http://www.t3db.ca/toxins/", chemIDs, ">", chemIDs, "</a>", sep = "")
  )
  curated_res$chemical_ID <- htmllink

  fname <- paste("Stage5_annotation_results", sep = "")
  unlink(fname)

  for (file in c("/stage2/", "/stage3/", "/stage4/", "/stage5/")) {
    outloc2 <- paste(outloc, file, sep = "")
    unlink(outloc2, force = TRUE, recursive = TRUE)
  }

  suppressWarnings(unlink("*.Rda"))

  for (file in c("step1_results.Rda", "plot.pdf", "Rplots.pdf", "Rplots.pdf")) {
    try(unlink(file), silent = TRUE)
  }

  curated_res <- as.data.frame(curated_res)
  curated_res <- curated_res[order(curated_res$Confidence, decreasing = TRUE), ]

  return(curated_res)
}
