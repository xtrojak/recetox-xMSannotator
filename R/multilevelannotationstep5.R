init.chemscoremat <- function(chemscoremat) {
  if (is.na(chemscoremat)) {
    chemscoremat <- read.csv("Stage4.csv")
  }
  chemscoremat <- as.data.frame(chemscoremat)
  chemscoremat$mz <- as.numeric(chemscoremat$mz)
  chemscoremat
}

multilevelannotationstep5 <- function(outloc,
                                      adduct_weights = NA,
                                      db_name = "HMDB",
                                      chemscoremat = NA,
                                      num_nodes = 2,
                                      scorethresh = 0) {
  setwd(outloc)
  curated_res <- init.chemscoremat(chemscoremat)





  if (is.na(adduct_weights)) {
    adduct_weights <- data.frame(Adduct = c("M+H", "M-H"),
                                 Weight = c(1, 1))
  }




  curated_res <- curated_res[order(curated_res$Confidence, curated_res$chemical_ID, curated_res$score, curated_res$Adduct, decreasing = TRUE), ]



  feature_count <- table(curated_res$mz)
  duplicated_features <- feature_count[which(feature_count > 1)]
  duplicated_features <- names(duplicated_features)

  bad_ind <- {}
  cl <- makeSOCKcluster(num_nodes)

  bad_ind <- foreach(mz_idx = 1:length(duplicated_features), .combine = rbind) %dopar% {
    mz <- duplicated_features[mz_idx]
    common_mz_idx <- which(curated_res$mz %in% mz)
    multimatch_features <- curated_res[common_mz_idx, ]

    features_adducts <- which(multimatch_features$Adduct %in% adduct_weights[, 1])
    if (length(features_adducts) > 0) {
      multimatch_features$score[features_adducts] <- (multimatch_features$score[features_adducts]) * 100
    }
    good_idx <- which(multimatch_features$score == max(multimatch_features$score, na.rm = TRUE))

    for (com_indval in 1:length(common_mz_idx)) {
      scoreval <- {}

      if (!com_indval %in% good_idx) {
        dmat_com <- curated_res[which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval]), ]
        scoreval <- ((dim(dmat_com)[1]) - 1) * dmat_com$score[1] / (dim(dmat_com)[1])
        scorevec <- c(rep(scoreval, length(which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval]))))

        if (length(scorevec) < length(which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval]))) {
          break
        }
        curated_res$score[which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[com_indval])] <- scorevec
      }
    }
    common_mz_idx <- common_mz_idx[-good_idx]

    return(common_mz_idx)
  }



  stopCluster(cl)

  if (length(bad_ind) > 0) {
    curated_res_unique <- curated_res[-c(bad_ind), ]
  } else {
    curated_res_unique <- curated_res
  }

  good_ind <- which(curated_res_unique$score >= scorethresh)

  curated_res_unique_highconf <- {}
  if (length(good_ind) > 0) {
    curated_res_unique_highconf <- curated_res_unique[good_ind, ]
  }

  curated_res <- curated_res_unique_highconf
  t2 <- table(curated_res$mz)

  same1 <- which(t2 == 1)

  uniquemz <- names(same1)

  curated_res$MatchCategory <- rep("Multiple", dim(curated_res)[1])

  curated_res$MatchCategory[which(curated_res$mz %in% uniquemz)] <- "Unique"


  write.csv(curated_res, file = "Stage5.csv", row.names = FALSE)

  curated_res$mz <- sprintf("%.5f", as.numeric(curated_res$mz))





  chemIDs <- curated_res$chemical_ID
  htmllink <- curated_res$chemical_ID

  link_text <- chemIDs[1]

  switch(db_name,
         "HMDB" = htmllink <- paste("<a href=http://www.hmdb.ca/metabolites/", chemIDs, ">", chemIDs, "</a>", sep = ""),
         "KEGG" = htmllink <- paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", chemIDs, ">", chemIDs, "</a>", sep = ""),
         "LipidMaps" = htmllink <- paste("<a href=http://www.lipidmaps.org/data/LMSDRecord.php?LMID=", chemIDs, ">", chemIDs, "</a>", sep = ""),
         "T3DB" = htmllink <- paste("<a href=http://www.t3db.ca/toxins/", chemIDs, ">", chemIDs, "</a>", sep = ""))

  fname <- paste("Stage5_annotation_results", sep = "")
  unlink(fname)

  curated_res$chemical_ID <- htmllink


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
