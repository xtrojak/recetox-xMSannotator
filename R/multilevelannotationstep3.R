compute_score <- function(chemscoremat, matrix, pathwaycheckmode, scorethresh, adduct_weights, max_diff_rt, bad_path_IDs, db_name) {
  pthresh = 0.05

  pathway_ids <- as.character(matrix[, 2])
  
  pathway_ids <- unique(pathway_ids)
  
  module_num <- gsub(chemscoremat$Module_RTclust,
                     pattern = "_[0-9]*",
                     replacement = "")
  
  chemscoremat <- cbind(chemscoremat, module_num)
  
  total_chem_count <- length(unique(matrix$chemid))
  
  if (! is.na(pathwaycheckmode)) {
    for (path_id in pathway_ids)
    {
      if (! path_id %in% bad_path_IDs) {
        pathway_chemicals <- matrix[which(matrix[, 2] %in% path_id), 1] 
        
        good_indices <- which(
          chemscoremat$chemical_ID %in% pathway_chemicals &
            chemscoremat$score >= scorethresh &
            chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
        )
        
        curmchemicaldata1 <- chemscoremat[good_indices, ]
        
        #molecules of interest in pathway (a)
        num_chems_inpath <- length(unique(curmchemicaldata1$chemical_ID))
        
        #total number of chemicals in pathway
        all_cur_path_numchem <- length(unique(pathway_chemicals))
        
        #non-focus molecules associated with pathway (c)
        num_chem_inpath_notinterest <- all_cur_path_numchem - num_chems_inpath
        
        indices <- which(
          chemscoremat$score >= scorethresh &
            chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
        )
        
        curmchemicaldata2 <- chemscoremat[indices, ]
        
        curmchemicaldata2 <- curmchemicaldata2[-which(curmchemicaldata2$chemical_ID %in% pathway_chemicals), ]
        
        #focus molecules not associated with this pathway (b)
        num_chems_notinpath <- length(unique(curmchemicaldata2$chemical_ID))
        
        all_notcurpath_numchem <- length(matrix[-which(matrix[, 1] %in% curmchemicaldata2$chemical_ID), 1])
        
        counts = matrix(
          data = c(
            num_chems_inpath,
            num_chem_inpath_notinterest,
            num_chems_notinpath,
            all_notcurpath_numchem
          ),
          nrow = 2
        )
        
        rm(curmchemicaldata2)
        p1 <- fisher.test(counts)
        p1 <- p1$p.value
        
        if (p1 > pthresh) {
          next
        } else {
          t1 <- table(curmchemicaldata1$module_num)
          
          module_counts <- t1[which(t1 > 0)]
          
          module_names <- names(module_counts)
          
          pathway_chemicals_1 <- curmchemicaldata1$chemical_ID
          
          # compatibility with wrong behavior
          if (db_name == "KEGG") {
            pathway_chemicals_to_iterate <- pathway_chemicals_1
          } else {
            pathway_chemicals_to_iterate <- pathway_chemicals
          }
          
          for (chemname in pathway_chemicals_to_iterate) {
            pathway_indices <- which(as.character(curmchemicaldata1$chemical_ID) == chemname)
            curmchemicaldata <- curmchemicaldata1[pathway_indices, ]
            chem_path_data <- matrix[which(matrix$chemid == chemname), ]
            t2 <- table(curmchemicaldata$module_num)
            cur_module <- names(t2[which(t2 == max(t2)[1])])
            
            mzid_cur <- paste(curmchemicaldata$mz,
                              curmchemicaldata$time,
                              sep = "_")
            
            if (nrow(curmchemicaldata) > 0) {
              pathwayscurchemical <- matrix[which(matrix[, 1] == chemname), 2]
              #for(cur_module in cur_modules)
              {
                if (pathwaycheckmode == "pm") {
                  num_chems <- t1[as.character(cur_module)]
                }
                
                num_chems <- round(num_chems, 0)
                
                module_indices <- which(curmchemicaldata$module_num == cur_module)
                num_chems_inmodule <- length(unique(curmchemicaldata1$chemical_ID[module_indices]))
                
                module_indices_in_pathway <- which(
                  chemscoremat$module_num == cur_module &
                    chemscoremat$chemical_ID %in% pathway_chemicals &
                    chemscoremat$score >= scorethresh &
                    chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                )
                cur_module_data <- chemscoremat[module_indices_in_pathway, ]
                a <- length(unique(cur_module_data$chemical_ID))
                rm(cur_module_data)
                
                module_indices <- which(
                  chemscoremat$module_num == cur_module &
                    chemscoremat$score >= scorethresh &
                    chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                )
                cur_module_data2 <- chemscoremat[module_indices, ]
                b <- length(unique(cur_module_data2$chemical_ID)) - a
                rm(cur_module_data2)
                
                other_module_indices_in_pathway <- which(
                  chemscoremat$module_num != cur_module &
                    chemscoremat$chemical_ID %in% pathway_chemicals &
                    chemscoremat$score >= scorethresh &
                    chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                )
                other_module_data <- chemscoremat[other_module_indices_in_pathway, ]
                c <- length(unique(other_module_data$chemical_ID))
                rm(other_module_data)
                
                other_module_indices <- which(
                  chemscoremat$module_num != cur_module &
                    chemscoremat$score >= scorethresh &
                    chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                )
                other_module_data2 <- chemscoremat[other_module_indices, ]
                d <- length(unique(other_module_data2$chemical_ID)) - c
                rm(other_module_data2)
                
                
                counts = matrix(data = c(a, c, b, d), nrow = 2)
                if (a > 1) {
                  p1 <- fisher.test(counts)
                  p1 <- p1$p.value
                } else {
                  p1 = 1
                }
                
                if (p1 > 0.2) {
                  next
                } else {
                  if (num_chems < 3) {
                    num_chems <- 0
                  } else {
                    if (is.na(curmchemicaldata$score[1])) {
                      diff_rt <- max(curmchemicaldata$time) - min(curmchemicaldata$time)
                      
                      if (diff_rt > max_diff_rt) {
                        if (length(which(t2 > 1)) == 1) {
                          curmchemicaldata$score <- rep(0.1, length(curmchemicaldata$score))
                        } else {
                          curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                        }
                      } else {
                        curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                      }
                    }
                    
                    # compatibility with wrong behavior
                    if (db_name == "KEGG") {
                      chemical_name <- c
                    } else {
                      chemical_name <- chemname
                    }
                    
                    if (curmchemicaldata$score[1] < scorethresh) {
                      indices <- which(
                        as.character(chemscoremat$chemical_ID) == chemical_name &
                          chemscoremat$Adduct %in% as.character(adduct_weights[, 1])
                      )
                      chemscoremat$score[indices] = as.numeric(chemscoremat$score[indices][1]) + num_chems
                    } else {
                      chemical_indices <- which(as.character(chemscoremat$chemical_ID) == chemical_name)
                      chemscoremat$score[chemical_indices] = chemscoremat$score[chemical_indices][1] + num_chems
                    }
                  }
                }
              }
            }
          }
        }
      }
      rm(curmchemicaldata1)
    }
  }
  
  chemscoremat <- replace_x_names(chemscoremat)
  return(chemscoremat)
}


compute_matrix <- function(DB, index, pattern) {
  temp_matrix <- apply(DB, 1, function(x) {
    chemid <- x[1]
    
    g1 <- gregexpr(x[index], pattern = pattern)
    regexp_check <- attr(g1[[1]], "match.length")
    if (regexp_check[1] < 0) {
      pathid = "-"
      
      return(cbind(chemid, pathid))
    } else {
      pathid <- strsplit(x = x[index], split = ";")
      
      pathid <- unlist(pathid)
      return(cbind(chemid, pathid))
    }
  })
  
  matrix <- ldply(temp_matrix, rbind)
  return(matrix)
}


load_chemscoremat <- function(num_sets) {
  chemscoremat <- lapply(1:num_sets, function(sind)
  {
    cur_fname <- paste("chem_score", sind, ".Rda", sep = "")
    load(cur_fname)
    
    curchemscoremat <- as.data.frame(curchemscoremat)
    return(curchemscoremat)
  })
  chemscoremat <- ldply(chemscoremat, rbind)
}


select_y_names <- function(chemscoremat, column_names) {
  #y because we want chemCompMZ ID and Name
  column_names_y <- column_names
  column_names_y[1] <- "cur_chem_score"
  column_names_y[7:8] <- paste(column_names_y[7:8], ".y", sep = "")
  column_names_y[11] <- "tempadduct"
  chemscoremat <- chemscoremat[, column_names_y]
  return(chemscoremat)
}


replace_x_names <- function(chemscoremat) {
  cnames <- colnames(chemscoremat)
  cnames <- gsub(cnames, pattern = ".x", replacement = "")
  colnames(chemscoremat) <- cnames
  return(chemscoremat)
}


sanitize_chemscoremat <- function(chemscoremat, chemCompMZ, column_names) {
  chemscoremat$Formula <- gsub(chemscoremat$Formula,
                               pattern = "_.*",
                               replacement = "")
  
  tempadduct <- chemscoremat$Adduct
  
  chemscoremat$Adduct <- gsub(chemscoremat$Adduct,
                              pattern = "_.*",
                              replacement = "")
  
  chemscoremat <- cbind(chemscoremat, tempadduct)
  
  chemscoremat <- merge(chemscoremat,
                        chemCompMZ[, c(2:4, 6)], 
                        by = c("Formula", "Adduct"))
  
  chemscoremat <- select_y_names(chemscoremat, column_names)
  colnames(chemscoremat) <- column_names
  return(chemscoremat)
}


multilevelannotationstep3 <- function(outloc1,
                                      chemscoremat = NA,
                                      adduct_weights = NA,
                                      num_sets = NA,
                                      db_name = NA,
                                      max_diff_rt = NA,
                                      pathwaycheckmode = "p",
                                      scorethresh = 0.1) {
  setwd(outloc1)
  
  load("chemCompMZ.Rda")

  outloc <- outloc1
  
  if (is.na(adduct_weights)) {
    adduct_weights <- data.frame(Adduct = c("M+H", "M-H"), Weight = c(1, 1))
  }
  
  outloc1 <- paste(outloc, "/stage2/", sep = "")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1)
  
  if (is.na(chemscoremat)) {
    chemscoremat <- load_chemscoremat(num_sets)
  }
  
  column_names <- c("score",
                    "Module_RTclust",
                    "mz",
                    "time",
                    "MatchCategory",
                    "theoretical.mz",
                    "chemical_ID",
                    "Name",
                    "Formula",
                    "MonoisotopicMass",
                    "Adduct",
                    "ISgroup",
                    "mean_int_vec",
                    "MD"
                   )
  
  chemscoremat <- sanitize_chemscoremat(chemscoremat, chemCompMZ, column_names)
  
  hmdbbad <- c("HMDB29244", "HMDB29245", "HMDB29246")
  bad_indices <- which(chemscoremat$chemical_ID %in% hmdbbad)
  
  if (length(bad_indices) > 0) {
    chemscoremat <- chemscoremat[-bad_indices, ]
  }
  
  if (db_name == "KEGG") {
    data(keggotherinf)
    
    matrix <- compute_matrix(keggotherinf, 4, "map")
    
    chemscoremat <- merge(chemscoremat,
                          keggotherinf,
                          by.x = "chemical_ID",
                          by.y = "KEGGID")
    
    chemscoremat <- compute_score(chemscoremat, matrix, pathwaycheckmode, scorethresh, adduct_weights, max_diff_rt, c("-", "map01100"), db_name)
  }
  else {
    if (db_name == "HMDB") {
      data(hmdbAllinf)
      
      hmdbAllinfv3.5 <- hmdbAllinf[, -c(26:27)]
      rm(hmdbAllinf, envir = .GlobalEnv)

      matrix <- compute_matrix(hmdbAllinfv3.5, 14, "SM")
      
      chemscoremat <- merge(chemscoremat,
                            hmdbAllinfv3.5,
                            by.x = "chemical_ID",
                            by.y = "HMDBID")
      
      rm(hmdbAllinfv3.5)
      
      chemscoremat <- compute_score(chemscoremat, matrix, pathwaycheckmode, scorethresh, adduct_weights, max_diff_rt, c("-"), db_name)
    }
  }
  
  chemscoremat <- replace_x_names(chemscoremat)

  good_ind <- which(chemscoremat$score >= scorethresh)
  if (length(good_ind) == 0) {
    chemscoremat <- {}
  } else {
    # rotate column names
    column_names[1:7] <- c(column_names[7], column_names[1:6])
    chemscoremat <- chemscoremat[, column_names]
  }
  
  write.csv(chemscoremat, file = "../Stage3.csv", row.names = FALSE)
  
  rm("outloc",
    "num_sets",
    "db_name",
    "num_sets",
    "adduct_weights",
    "chemCompMZ",
    "hmdbAllinf",
    "hmdbAllinfv3.6")
  
  return(chemscoremat)
}
