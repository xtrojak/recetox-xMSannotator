multilevelannotationstep5 <- function(outloc,
                                      adduct_weights = NA,
                                      db_name = "HMDB",
                                      chemscoremat = NA,
                                      num_nodes = 2) {
  setwd(outloc)

  if (is.na(chemscoremat)) {
    curated_res <- read.csv("Stage4.csv")
  } else {
    curated_res <- chemscoremat
    rm(chemscoremat)
  }
  curated_res <- as.data.frame(curated_res)
  curated_res$mz <- as.numeric(curated_res$mz)

  scorethresh <- 0




  if (is.na(adduct_weights)) {
    adduct_weights <- data.frame(Adduct = c("M+H", "M-H"),
                                 Weight = c(1, 1))
  }




  curated_res <- curated_res[order(curated_res$Confidence, curated_res$chemical_ID, curated_res$score, curated_res$Adduct, decreasing = TRUE), ]




  t1 <- table(curated_res$mz, curated_res$chemical_ID)

  s1 <- apply(t1, 1, sum)

  rm(t1)


  s2 <- s1[which(s1 > 1)]
  mzdup <- names(s2)
  mzdup <- as.data.frame(mzdup)
  mzdup <- mzdup[, 1]

  bad_ind <- {}
  cl <- makeSOCKcluster(num_nodes)

  bad_ind <- foreach(mind = 1:length(mzdup), .combine = rbind) %dopar% {
    mznum <- mzdup[mind]
    dmultsub <- curated_res[which(curated_res$mz %in% mznum), ]
    dgood_add <- which(dmultsub$Adduct %in% adduct_weights[, 1])
    if (length(dgood_add) > 0) {
      dmultsub$score[dgood_add] <- (dmultsub$score[dgood_add]) * 100
    }
    com_ind <- which(curated_res$mz %in% mznum)
    good_ind <- which(dmultsub$score == max(dmultsub$score, na.rm = TRUE))

    for (com_indval in 1:length(com_ind)) {
      scoreval <- {}
      if (com_indval %in% good_ind == FALSE) {
        dmat_com <- curated_res[which(curated_res$chemical_ID %in% dmultsub$chemical_ID[com_indval]), ]
        scoreval <- ((dim(dmat_com)[1]) - 1) * dmat_com$score[1] / (dim(dmat_com)[1])
        scorevec <- c(rep(scoreval, length(which(curated_res$chemical_ID %in% dmultsub$chemical_ID[com_indval]))))
        if (length(scorevec) < length(which(curated_res$chemical_ID %in% dmultsub$chemical_ID[com_indval]))) {
          break
        }
        curated_res$score[which(curated_res$chemical_ID %in% dmultsub$chemical_ID[com_indval])] <- scorevec
      }
    }
    com_ind <- com_ind[-good_ind]

    return(com_ind)
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

  if (db_name == "HMDB") {
    htmllink <- paste("<a href=http://www.hmdb.ca/metabolites/", chemIDs, ">", chemIDs, "</a>", sep = "")
  } else {
    if (db_name == "KEGG") {
      htmllink <- paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", chemIDs, ">", chemIDs, "</a>", sep = "")
    } else {
      if (db_name == "LipidMaps") {
        htmllink <- paste("<a href=http://www.lipidmaps.org/data/LMSDRecord.php?LMID=", chemIDs, ">", chemIDs, "</a>", sep = "")
      } else {
        if (db_name == "T3DB") {
          htmllink <- paste("<a href=http://www.t3db.ca/toxins/", chemIDs, ">", chemIDs, "</a>", sep = "")
        }
      }
    }
  }

  fname <- paste("Stage5_annotation_results", sep = "")
  unlink(fname)

  curated_res$chemical_ID <- htmllink



  outloc2 <- paste(outloc, "/stage2/", sep = "")

  unlink(outloc2, force = TRUE, recursive = TRUE)

  outloc2 <- paste(outloc, "/stage2", sep = "")

  unlink(outloc2, force = TRUE, recursive = TRUE)

  file.remove(dir(outloc2, full.names = TRUE))

  outloc2 <- paste(outloc, "/stage3/", sep = "")
  unlink(outloc2, force = TRUE, recursive = TRUE)
  file.remove(dir(outloc2, full.names = TRUE))
  outloc2 <- paste(outloc, "/stage4/", sep = "")
  file.remove(dir(outloc2, full.names = TRUE))
  unlink(outloc2, force = TRUE, recursive = TRUE)
  outloc2 <- paste(outloc, "/stage5/", sep = "")
  file.remove(dir(outloc2, full.names = TRUE))
  unlink(outloc2, force = TRUE, recursive = TRUE)

  suppressWarnings(unlink("*.Rda"))


  try(unlink("step1_results.Rda"), silent = TRUE)
  try(unlink("plot.pdf"), silent = TRUE)
  try(unlink("Rplots.pdf"), silent = TRUE)
  try(unlink("Rplots.pdf"), silent = TRUE)

  curated_res$chemical_ID <- chemIDs

  curated_res <- as.data.frame(curated_res)

  curated_res <- curated_res[order(curated_res$Confidence, decreasing = TRUE), ]

  print("Stage 5 confidence level distribution for unique chemical/metabolite IDs")
  print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID) == TRUE)]))


  print("Stage 5 confidence level distribution for unique chemical/metabolite formulas")
  print(table(curated_res$Confidence[-which(duplicated(curated_res$Formula) == TRUE)]))


  return(curated_res)
}
