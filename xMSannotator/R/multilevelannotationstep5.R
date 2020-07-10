multilevelannotationstep5 <- function(outloc, max.mz.diff = 5, max.rt.diff = 30, adduct_table, adduct_weights, filter.by = NA, min_ions_perchem = 1, boostIDs = NA, max_isp = 5, dbAllinf = NA, db_name = "HMDB", chemscoremat = NA, num_nodes = 2) {
  setwd(outloc)

  if (is.na(chemscoremat) == TRUE) {
    curated_res <- read.csv("Stage4.csv")
  }else {
    curated_res <- chemscoremat
    rm(chemscoremat)
  }
  curated_res <- as.data.frame(curated_res)
  curated_res$mz <- as.numeric(curated_res$mz)

  scorethresh <- 0
  adduct_table <- adduct_table[order(adduct_table$Adduct),]

  t1 <- table(curated_res$mz, curated_res$chemical_ID)
  t2 <- apply(t1, 1, sum)

  multi_mz <- names(t2[which(t2 > 1)])

  curated_res$MatchCategory <- gsub(as.character(curated_res$MatchCategory), pattern = "Multiple", replacement = "Unique")
  curated_res$MatchCategory[which(curated_res$mz %in% multi_mz)] <- "Multiple"

  curated_res <- curated_res[order(curated_res$Confidence, curated_res$chemical_ID, curated_res$score, curated_res$Adduct, decreasing = TRUE),]

  t1 <- table(curated_res$mz, curated_res$chemical_ID)
  s1 <- apply(t1, 1, sum)

  rm(t1)

  mzunique <- colnames(which(s1 == 1))

  dunique <- curated_res[which(curated_res$mz %in% mzunique),]

  s2 <- s1[which(s1 > 1)]
  mzdup <- names(s2)
  mzdup <- as.data.frame(mzdup)
  mzdup <- mzdup[, 1]

  cl <- parallel::makeSOCKcluster(num_nodes)

  bad_ind <- foreach::foreach(mind = seq_along(mzdup), .combine = rbind) %dopar% {
    mznum <- mzdup[mind]
    dmultsub <- curated_res[which(curated_res$mz %in% mznum),]
    dgood_add <- which(dmultsub$Adduct %in% adduct_weights[, 1])
    if (length(dgood_add) > 0) {
      dmultsub$score[dgood_add] <- (dmultsub$score[dgood_add]) * 100
    }
    com_ind <- which(curated_res$mz %in% mznum)
    good_ind <- which(dmultsub$score == max(dmultsub$score, na.rm = TRUE))

    for (com_indval in seq_along(com_ind)) {
      if (com_indval %in% good_ind == FALSE) {
        bla <- which(curated_res$chemical_ID %in% dmultsub$chemical_ID[com_indval])
        dmat_com <- curated_res[bla,]
        scoreval <- ((dim(dmat_com)[1]) - 1) * dmat_com$score[1] / (dim(dmat_com)[1])
        curated_res$score[bla] <- rep(scoreval, length(bla))
      }
    }
    com_ind <- com_ind[-good_ind]
    return(com_ind)
  }

  parallel::stopCluster(cl)

  curated_res <- curated_res[-bad_ind,]
  curated_res <- curated_res[curated_res$score >= scorethresh,]
  curated_res$MatchCategory <- factor(duplicated(curated_res$mz), labels = c('Unique', 'Multiple'))

  write.csv(curated_res, file = "Stage5.csv", row.names = FALSE)

  curated_res$mz <- sprintf("%.5f", as.numeric(curated_res$mz))

  curated_res <- as.data.frame(curated_res)
  curated_res <- curated_res[order(curated_res$Confidence, decreasing = TRUE),]

  print("Stage 5 confidence level distribution for unique chemical/metabolite IDs")
  print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID) == TRUE)]))

  print("Stage 5 confidence level distribution for unique chemical/metabolite formulas")
  print(table(curated_res$Confidence[-which(duplicated(curated_res$Formula) == TRUE)]))

  return(curated_res)
}