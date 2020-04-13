multilevelannotation <- function(dataA, max.mz.diff = 10, max.rt.diff = 10, cormethod = "pearson", num_nodes = parallel::detectCores(), queryadductlist = c("all"), gradienttype = "Acetonitrile", mode = "pos", outloc, db_name = "HMDB", adduct_weights = NA, num_sets = 3000, allsteps = TRUE,
    corthresh = 0.7, NOPS_check = TRUE, customIDs = NA, missing.value = NA, deepsplit = 2, networktype = "unsigned", minclustsize = 10, filter.by = c("M+H"), redundancy_check = TRUE, min_ions_perchem = 1, biofluid.location = NA, origin = NA,
    status = NA, boostIDs = NA, max_isp = 5, customDB = NA, HMDBselect = "union", pathwaycheckmode = "pm") {
    
    options(warn = -1)
    
    WGCNA::enableWGCNAThreads(nThreads = num_nodes)
    
    dataA <- as.data.frame(dataA)
    dataA$mz <- round(dataA$mz, 5)
    dataA$time <- round(dataA$time, 1)
    
    dataA[, -(1:2)] <- round(dataA[, -(1:2)], 1)
    
    if (is.na(customIDs) == FALSE) {
        customIDs <- as.data.frame(customIDs)
    }
    
    data(adduct_table)
    
    print(paste0("Annotating using ", db_name, " database:"))
    
    max_diff_rt <- max.rt.diff
    cutheight <- 1 - corthresh  #module.merge.dissimilarity
    time_step <- 1
    step1log2scale <- FALSE
    
    adduct_table <- unique(adduct_table)
    
    suppressWarnings(dir.create(outloc))

    current_dir <- getwd()
    on.exit(setwd(current_dir))
    setwd(outloc)
    
    if (queryadductlist == "all" & mode == "pos") {
        adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "positive") | (adduct_table$Type == gradienttype & adduct_table$Mode == "positive")]
        adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
    } else {
        if (queryadductlist == "all" & mode == "neg") {
            adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "negative") | (adduct_table$Type == gradienttype & adduct_table$Mode == "negative")]
            adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
        } else {
            adduct_names <- adduct_table$Adduct[which(adduct_table$Adduct %in% queryadductlist)]
            adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
        }
    }
    
    outloc_allres <- outloc
    
    l1 <- list.files(outloc_allres)
    check_step1 <- which(l1 == "step1_results.Rda")
    
    if (length(check_step1) < 1) {
        print("Status 1: Computing modules using WGCNA")
        
        check_levelA <- which(l1 == "xMSannotator_levelA_modules.Rda")
        
        if (is.na(missing.value) == FALSE) {
            dataA <- replace(as.matrix(dataA), which(dataA == missing.value), NA)
            dataA <- as.data.frame(dataA)
        }
        
        dataA <- unique(dataA)
        mzid <- paste(dataA$mz, dataA$time, sep = "_")
        
        if (length(which(duplicated(mzid) == TRUE)) > 0) {
            dataA <- dataA[-which(duplicated(mzid) == TRUE), ]
        }
        
        system.time(global_cor <- WGCNA::cor(t(dataA[, -(1:2)]), nThreads = num_nodes, method = cormethod, use = "p"))
        
        global_cor <- round(global_cor, 2)
        
        dataA <- unique(dataA)
        setwd(outloc_allres)
        
        mzid <- paste(dataA$mz, dataA$time, sep = "_")
        colnames(global_cor) <- mzid
        rownames(global_cor) <- mzid
        
        save(global_cor, file = "global_cor.Rda")
        
        hr <- flashClust::flashClust(as.dist(1 - global_cor), method = "complete")
        dissTOMCormat <- (1 - global_cor)
        
        a1 <- apply(global_cor, 1, function(x) {
            length(which(x > corthresh))
        })
        
        fname_c <- paste0("NumberOfCorrelationsPerFeature_cor", corthresh, ".csv")
        
        write.csv(a1, file = fname_c)
        
        rm(global_cor)
        
        # using dynamic
        mycl_metabs <- dynamicTreeCut::cutreeDynamic(hr, distM = dissTOMCormat, deepSplit = 1, minClusterSize = minclustsize, pamRespectsDendro = FALSE, pamStage = TRUE, verbose = 0)
        
        clustmethod <- "WGCNA"
        if (length(check_levelA) < 1) {
            if (clustmethod == "WGCNA") {
                levelA_res <- get_peak_blocks_modulesvhclust(dataA = dataA, simmat = NA, adjacencyfromsimilarity = FALSE, time_step = time_step, max.rt.diff = max_diff_rt, outloc, column.rm.index = NA, cor.thresh = NA, deepsplit = deepsplit, minclustsize = minclustsize, 
                  cutheight = cutheight, cormethod = cormethod, networktype = networktype, num_nodes = num_nodes, step1log2scale = step1log2scale, mycl_metabs = mycl_metabs)
                
                setwd(outloc_allres)
                levelA_res <- levelA_res[, -(1:4)]
                save(levelA_res, file = "xMSannotator_levelA_modules.Rda")
            } else {
                if (clustmethod == "graph") {
                  levelA_res <- get_peak_blocks_graph(dataA, simmat = global_cor, adjacencyfromsimilarity = TRUE, time_step = 3, max.rt.diff = max_diff_rt, outloc, column.rm.index = NA, cor.thresh = NA, cormethod = cormethod, networktype = networktype, num_nodes = num_nodes)
                } else {
                  if (clustmethod == "hclust") {
                    levelA_res <- get_peak_blocks_hclust(dataA, time_step = time_step, max.rt.diff = max_diff_rt, outloc, column.rm.index = NA, cor.thresh = NA, deepsplit = deepsplit, minclustsize = minclustsize, cutheight = cutheight, cormethod = cormethod)
                  }
                }
            }
            
            setwd(outloc_allres)
            write.csv(levelA_res, file = "Stage1.csv", row.names = FALSE)
        } else {
            setwd(outloc_allres)
            load("xMSannotator_levelA_modules.Rda")
        }
        
        setwd(outloc)
        
        levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ]
        levelA_res <- levelA_res[, 1:3]
        
        dataA <- dataA[order(dataA$mz, dataA$time), ]
        
        mean_int_vec <- apply(dataA[, -(1:2)], 1, function(x) {
            mean(x, na.rm = TRUE)
        })
        
        dataA <- dataA[, 1:2]
        
        if (queryadductlist == "all" & mode == "pos") {
            adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "positive") | (adduct_table$Type == gradienttype & adduct_table$Mode == "positive")]
            adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
        } else {
            if (queryadductlist == "all" & mode == "neg") {
                adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "negative") | (adduct_table$Type == gradienttype & adduct_table$Mode == "negative")]
                adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
            } else {
                adduct_names <- adduct_table$Adduct[which(adduct_table$Adduct %in% queryadductlist)]
                adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
            }
        }
        
        adduct_names <- unique(adduct_names)
        
        if (db_name == "HMDB") {
            data(hmdbAllinf)
            
            hmdbAllinf$Name <- gsub(hmdbAllinf$Name, pattern = "[\\\"']", replacement = "")
            hmdbAllinfv3.6 <- hmdbAllinf
            
            rm(hmdbAllinf)
            rm(hmdbAllinf, envir = .GlobalEnv)
            
            suppressWarnings(if (is.na(customIDs) == TRUE) {
                if (is.na(biofluid.location) == FALSE & is.na(origin) == TRUE) {
                  gres <- gregexpr(hmdbAllinfv3.6$BioFluidLocation, pattern = biofluid.location, ignore.case = TRUE)
                  
                  if (is.na(status) == FALSE) {
                    gres3 <- gregexpr(hmdbAllinfv3.6$HMDBStatus, pattern = status, ignore.case = TRUE)
                    
                    gres3_good <- which(gres3 == 1)
                    gres_good <- which(gres == 1)
                    if (HMDBselect == "intersect") {
                      gres <- intersect(gres_good, gres3_good)
                    } else {
                      gres <- c(gres_good, gres3_good)
                      gres <- unique(gres)
                    }
                  } else {
                    gres <- which(gres == 1)
                  }
                  
                  customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]
                  rm(hmdbAllinfv3.6)
                } else {
                  if (is.na(biofluid.location) == FALSE & is.na(origin) == FALSE) {
                    gres <- gregexpr(hmdbAllinfv3.6$BioFluidLocation, pattern = biofluid.location, ignore.case = TRUE)
                    gres2 <- gregexpr(hmdbAllinfv3.6$Origin, pattern = origin, ignore.case = TRUE)
                    
                    gres_good <- which(gres == 1)
                    gres2_good <- which(gres2 == 1)
                    
                    if (is.na(status) == FALSE) {
                      gres3 <- gregexpr(hmdbAllinfv3.6$HMDBStatus, pattern = status, ignore.case = TRUE)
                      gres3_good <- which(gres3 == 1)
                      
                      if (HMDBselect == "intersect") {
                        gres <- intersect(gres_good, gres2_good, gres3_good)
                      } else {
                        gres <- c(gres_good, gres2_good, gres3_good)
                        gres <- unique(gres)
                      }
                    } else {
                      if (HMDBselect == "intersect") {
                        gres <- intersect(gres_good, gres2_good)
                      } else {
                        gres <- c(gres_good, gres2_good)
                        gres <- unique(gres)
                      }
                    }
                    
                    customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]
                  } else {
                    if (is.na(biofluid.location) == TRUE & is.na(origin) == FALSE) {
                      gres <- gregexpr(hmdbAllinfv3.6$Origin, pattern = origin, ignore.case = TRUE)
                      
                      if (is.na(status) == FALSE) {
                        gres3 <- gregexpr(hmdbAllinfv3.6$HMDBStatus, pattern = status, ignore.case = TRUE)
                        gres3_good <- which(gres3 == 1)
                        gres_good <- which(gres == 1)
                        
                        if (HMDBselect == "intersect") {
                          gres <- intersect(gres_good, gres3_good)
                        } else {
                          gres <- c(gres_good, gres3_good)
                          gres <- unique(gres)
                        }
                      } else {
                        gres <- which(gres == 1)
                      }
                      
                      customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]
                    } else {
                      if (is.na(status) == FALSE) {
                        gres3 <- gregexpr(hmdbAllinfv3.6$HMDBStatus, pattern = status, ignore.case = TRUE)
                        gres3_good <- which(gres3 == 1)
                        gres <- gres3_good
                        
                        customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]
                        
                      } else {
                        customIDs <- hmdbAllinfv3.6[, c(1, 20)]
                      }
                    }
                  }
                }
            })
            
            rm(hmdbAllinfv3.6, envir = .GlobalEnv)
            data(hmdbCompMZ)
            
            hmdbCompMZ$mz <- round(as.numeric(as.character(hmdbCompMZ$mz)), 5)
            hmdbCompMZ$Name <- gsub(hmdbCompMZ$Name, pattern = "[\\\"']", replacement = "")
            
            suppressWarnings(if (is.na(customIDs) == FALSE) {
                customIDs <- unique(customIDs)
                hmdbCompMZ <- hmdbCompMZ[which(hmdbCompMZ$HMDBID %in% customIDs[, 1]), ]
                hmdbCompMZ <- hmdbCompMZ[which(hmdbCompMZ$Adduct %in% adduct_names), ]
            })
            
            hmdbCompMZ <- hmdbCompMZ[which(hmdbCompMZ$Adduct %in% adduct_names), ]
            chemCompMZ <- hmdbCompMZ
            
            print("Dimension of precomputed HMDB m/z database")
            print(dim(chemCompMZ))
            
            hmdbCompMZ <- 1
            hmdbAllinf <- 1
            try(rm(hmdbCompMZ), silent = TRUE)
            try(rm(hmdbCompMZ, envir = .GlobalEnv), silent = TRUE)
            try(rm(hmdbAllinf, envir = .GlobalEnv), silent = TRUE)
        } else {
            if (db_name == "KEGG") {
                data(keggCompMZ)
                
                keggCompMZ$mz <- round(as.numeric(as.character(keggCompMZ$mz)), 5)
                keggCompMZ$Name <- gsub(keggCompMZ$Name, pattern = "[\\\"']", replacement = "")
                
                suppressWarnings(if (is.na(customIDs) == FALSE) {
                  customIDs <- unique(customIDs)
                  keggCompMZ <- keggCompMZ[which(keggCompMZ$KEGGID %in% customIDs[, 1]), ]
                  keggCompMZ <- keggCompMZ[which(keggCompMZ$Adduct %in% adduct_names), ]
                })
                
                keggCompMZ <- keggCompMZ[which(keggCompMZ$Adduct %in% adduct_names), ]
                chemCompMZ <- keggCompMZ
                
                print("Dimension of precomputed KEGG m/z database")
                print(dim(chemCompMZ))
                
                rm(keggCompMZ)
                rm(keggCompMZ, envir = .GlobalEnv)
            } else {
                if (db_name == "LipidMaps") {
                  data(lipidmapsCompMZ)
                  
                  lipidmapsCompMZ <- lipidmapsCompMZ[which(lipidmapsCompMZ$Adduct %in% adduct_names), ]
                  lipidmapsCompMZ$mz <- round(as.numeric(as.character(lipidmapsCompMZ$mz)), 5)
                  lipidmapsCompMZ$Name <- gsub(lipidmapsCompMZ$Name, pattern = "[\\\"']", replacement = "")
                  chemCompMZ <- lipidmapsCompMZ
                  
                  print("Dimension of precomputed LipidMaps m/z database")
                  print(dim(chemCompMZ))
                  
                  rm(lipidmapsCompMZ)
                  rm(lipidmapsCompMZ, envir = .GlobalEnv)
                } else {
                  if (db_name == "T3DB") {
                    data(t3dbCompMZ)
                    
                    t3dbCompMZ <- t3dbCompMZ[which(t3dbCompMZ$Adduct %in% adduct_names), ]
                    t3dbCompMZ$mz <- round(as.numeric(as.character(t3dbCompMZ$mz)), 5)
                    t3dbCompMZ$Name <- gsub(t3dbCompMZ$Name, pattern = "[\\\"']", replacement = "")
                    chemCompMZ <- t3dbCompMZ
                    
                    print("Dimension of precomputed T3DB m/z database")
                    print(dim(chemCompMZ))
                    
                    rm(t3dbCompMZ)
                    rm(t3dbCompMZ, envir = .GlobalEnv)
                  } else {
                    if (db_name == "Custom") {
                      data(adduct_table)
                      
                      inputmassmat <- customDB
                      
                      if (length(which(duplicated(inputmassmat[, 1]) == TRUE)) > 0) {
                        inputmassmat <- inputmassmat[-which(duplicated(inputmassmat[, 1]) == TRUE), ]
                      }
                      
                      mz_search_list <- lapply(1:dim(inputmassmat)[1], function(m) {
                        adduct_names <- as.character(adduct_names)
                        
                        mz_search_list <- xMSannotator::get_mz_by_monoisotopicmass(monoisotopicmass = as.numeric(as.character(inputmassmat[m, 4])), dbid = inputmassmat[m, 1], name = as.character(inputmassmat[m, 2]), formula = as.character(inputmassmat[m, 3]), queryadductlist = adduct_names,
                          adduct_table = adduct_table)
                        return(mz_search_list)
                      })
                      
                      mz_search_list <- plyr::ldply(mz_search_list, rbind)
                      save(mz_search_list, file = "mz_search_list.Rda")
                      chemCompMZ <- mz_search_list
                      rm(mz_search_list)
                      rm(customDB)
                      rm(customDB, envir = .GlobalEnv)
                    } else {
                      stop("db_name should be: KEGG, HMDB, T3DB, LipidMaps, or Custom")
                    }
                  }
                }
            }
        }
        data(adduct_table)
        
        if (is.na(adduct_weights) == TRUE) {
            data(adduct_weights)
            adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
            adduct_weights1[1, ] <- c("M+H", 1)
            adduct_weights1[2, ] <- c("M-H", 1)
            adduct_weights <- as.data.frame(adduct_weights1)
            colnames(adduct_weights) <- c("Adduct", "Weight")
            adduct_weights <- as.data.frame(adduct_weights)
        }
        
        chemCompMZ$Name <- gsub("([-:();])|[[:punct:]]", "\\1", chemCompMZ$Name)
        chemCompMZ$Formula <- gsub("([-:();])|[[:punct:]]", "\\1", chemCompMZ$Formula)
        
        cnames <- colnames(chemCompMZ)
        cnames[2] <- "chemical_ID"
        colnames(chemCompMZ) <- cnames
        
        chemCompMZ <- as.data.frame(chemCompMZ)
        
        setwd(outloc)
        l1 <- list.files(outloc)
        check_levelB <- which(l1 == "xMSannotator_levelB.Rda")
        
        save(chemCompMZ, file = "chemCompMZ.Rda")
        
        if (length(check_levelB) < 1) {
            cl <- parallel::makeCluster(num_nodes)
            parallel::clusterEvalQ(cl, library(XML))
            parallel::clusterEvalQ(cl, library(R2HTML))
            parallel::clusterEvalQ(cl, library(RCurl))
            parallel::clusterEvalQ(cl, library(SSOAP))
            parallel::clusterEvalQ(cl, library(limma))
            parallel::clusterEvalQ(cl, library(plyr))
            parallel::clusterEvalQ(cl, "processWSDL")
            parallel::clusterEvalQ(cl, library(png))
            parallel::clusterExport(cl, "Annotationbychemical_IDschild_multilevel")
            parallel::clusterExport(cl, "Annotationbychemical_IDschild")
            parallel::clusterExport(cl, "find.Overlapping.mzs")
            parallel::clusterExport(cl, "find.Overlapping.mzsvparallel")
            parallel::clusterExport(cl, "overlapmzchild")
            parallel::clusterExport(cl, "getVenn")
            
            if (length(which(duplicated(chemCompMZ$Formula) == TRUE)) > 0) {
                if (db_name == "Custom") {
                  chemCompMZ$mz <- as.numeric(as.character(chemCompMZ$mz))
                }
                chemCompMZ_unique_formulas <- chemCompMZ[-which(duplicated(chemCompMZ$Formula) == TRUE), ]
                chemCompMZ_unique_formulas <- rbind(chemCompMZ_unique_formulas, chemCompMZ[which(chemCompMZ$chemical_ID %in% chemCompMZ_unique_formulas$chemical_ID), ])
                chemCompMZ_unique_formulas <- unique(chemCompMZ_unique_formulas)
                chemCompMZ_unique_formulas <- chemCompMZ_unique_formulas[order(chemCompMZ_unique_formulas$chemical_ID), ]
            } else {
                chemCompMZ_unique_formulas <- chemCompMZ
            }
            
            save(chemCompMZ, file = "chemCompMZ.Rda")
            rm(chemCompMZ)
            
            formula_table <- table(chemCompMZ_unique_formulas$Formula)
            uniq_formulas <- names(formula_table)
            formula_ID <- paste("Formula", seq(seq_along(uniq_formulas)), sep = "_")
            
            formula_id_mat <- cbind(formula_ID, uniq_formulas)
            formula_id_mat <- as.data.frame(formula_id_mat)
            colnames(formula_id_mat) <- c("Formula_ID", "Formula")
            
            chemCompMZ_unique_formulas <- merge(chemCompMZ_unique_formulas, formula_id_mat, by = "Formula")
            chemCompMZ_unique_formulas$chemical_ID <- chemCompMZ_unique_formulas$Formula_ID
            chemCompMZ_unique_formulas <- chemCompMZ_unique_formulas[, c("mz", "chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct", "AdductMass")]
            chemCompMZ_unique_formulas <- as.data.frame(chemCompMZ_unique_formulas)
            
            s1 <- seq(1, length(adduct_names))
            print("Status 2: Mapping m/z to metabolites:")
            
            adduct_names <- as.character(adduct_names)
            
            l2 <- parallel::parLapply(cl, s1, Annotationbychemical_IDschild, dataA = dataA, queryadductlist = adduct_names, adduct_type = c("S", gradienttype), adduct_table = adduct_table, max.mz.diff = max.mz.diff, outloc = outloc, keggCompMZ = chemCompMZ_unique_formulas,
                otherdbs = FALSE, otherinfo = FALSE, num_nodes = 1)
            
            parallel::stopCluster(cl)
            
            rm(chemCompMZ)
            levelB_res <- {
            }
            for (j in seq_along(l2)) {
                if (length(l2[[j]]) > 1) {
                  levelB_res <- rbind(levelB_res, l2[[j]])
                }
            }
            
            rm(l2)
            
            if (nrow(levelB_res) < 1) {
                stop("No matches found.")
            }
            
            MatchCategory <- rep("Multiple", dim(levelB_res)[1])
            
            levelB_res$mz <- as.numeric(as.character(levelB_res$mz))
            levelB_res$time <- as.numeric(as.character(levelB_res$time))
            levelB_res <- as.data.frame(levelB_res)
            levelB_res <- cbind(levelB_res, MatchCategory)
            
            print("DB matches")
            levelB_res <- unique(levelB_res)
            print(dim(levelB_res))
            
            uniq_formula <- as.character(unique(levelB_res$Formula))
            bad_formula <- which(is.na(uniq_formula) == TRUE)
            if (length(bad_formula) > 0) {
                uniq_formula <- uniq_formula[-bad_formula]
            }
            
            cl <- parallel::makeCluster(num_nodes)
            
            parallel::clusterExport(cl, "check_golden_rules")
            parallel::clusterExport(cl, "check_element")
            
            levelB_res_check <- parallel::parLapply(cl, seq_along(uniq_formula), function(j, uniq_formula, NOPS_check) {
                curformula <- as.character(uniq_formula[j])
                return(xMSannotator::check_golden_rules(curformula, NOPS_check = NOPS_check))
            }, uniq_formula = uniq_formula, NOPS_check = NOPS_check)
            parallel::stopCluster(cl)
            
            
            levelB_res_check2 <- plyr::ldply(levelB_res_check, rbind)
            levelB_res_check3 <- levelB_res_check2[which(levelB_res_check2[, 2] == 1), ]
            levelB_res <- levelB_res[which(levelB_res$Formula %in% as.character(levelB_res_check3[, 1])), ]
            
            water_adducts <- c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
            water_adduct_ind <- which(levelB_res$Adduct %in% water_adducts)
            
            cl <- parallel::makeCluster(num_nodes)
            parallel::clusterExport(cl, "check_element")
            
            if (length(water_adduct_ind) > 0) {
                levelB_res2 <- levelB_res[water_adduct_ind, ]
                levelB_res <- levelB_res[-water_adduct_ind, ]
                sind1 <- seq(1:dim(levelB_res2)[1])
                levelB_res_check3 <- parallel::parLapply(cl, sind1, function(j) {
                  curformula <- as.character(levelB_res2$Formula[j])
                  numoxygens <- xMSannotator::check_element(curformula, "O")
                  
                  if (numoxygens > 0) {
                    bool_check <- 1
                  } else {
                    bool_check <- 0
                  }
                  
                  res <- cbind(curformula, bool_check)
                  res <- as.data.frame(res)
                  return(res)
                })
                
                levelB_res_check4 <- plyr::ldply(levelB_res_check3, rbind)
                valid_form <- {
                }
                
                if (length(which(levelB_res_check4[, 2] == 1)) > 0) {
                  levelB_res_check4 <- levelB_res_check4[which(levelB_res_check4[, 2] == 1), ]
                  valid_form <- which(levelB_res2$Formula %in% as.character(levelB_res_check4[, 1]))
                }
                if (length(valid_form) > 0) {
                  levelB_res2 <- levelB_res2[valid_form, ]
                  levelB_res <- rbind(levelB_res, levelB_res2)
                }
            }
            
            colnames(levelB_res) <- c("theoretical.mz", "chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct", "AdductMass", "mz", "time", "MatchCategory")
            save(levelB_res, file = "xMSannotator_levelB.Rda")
        } else {
            print("Status 2: Using existing m/z mapping results:")
            load("xMSannotator_levelB.Rda")
        }
        
        try(rm(hmdbAllinf, envir = .GlobalEnv), silent = TRUE)
        try(rm(hmdbAllinfv3.6), silent = TRUE)
        try(rm(hmdbCompMZ), silent = TRUE)
        
        levelA_res$mz <- round(levelA_res$mz, 5)
        levelB_res$mz <- round(levelB_res$mz, 5)
        
        levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ]
        
        if (length(which(duplicated(mzid) == TRUE)) > 0) {
            levelA_res <- levelA_res[-which(duplicated(mzid) == TRUE), ]
            dataA <- dataA[-which(duplicated(mzid) == TRUE), ]
        }
        
        levelA_res <- levelA_res[order(levelA_res$Module_RTclust), ]
        module_num <- gsub(levelA_res$Module_RTclust, pattern = "_[0-9]{1,}", replacement = "")
        levelA_res_all <- levelA_res[, 1:2]
        levelA_res_all$Module_RTclust <- module_num
        mzdefect <- 1 * ((levelA_res$mz - floor(levelA_res$mz)))
        
        levelA_res2 <- cbind(mzdefect, levelA_res)
        
        massdefect_cor_groups <- sapply(list(myData1 = levelA_res2), function(x) split(x, cut(levelA_res2$mzdefect, breaks = seq(0, 1, 0.01))))
        
        if (length(which(levelA_res2$mzdefect == 0) > 0)) {
            massdefect_cor_groups[[length(massdefect_cor_groups) + 1]] <- levelA_res2[which(levelA_res2$mzdefect == 0), ]
        }
        diffmatB <- lapply(seq_along(massdefect_cor_groups), function(gnum) {
            cur_group <- {
            }
            
            if (dim(massdefect_cor_groups[[gnum]])[1] > 0) {
                ISgroup <- paste("ISgroup", massdefect_cor_groups[[gnum]]$Module_RTclust, gnum, sep = "_")
                
                cur_group <- as.data.frame(massdefect_cor_groups[[gnum]])
                cur_group <- cbind(ISgroup, cur_group)
                
                if (length(cur_group) > 0) {
                  cur_group <- cur_group[order(cur_group$mz, cur_group$time), ]
                }
            } else {
                print(gnum)
            }
            return(cur_group)
        })
        
        diffmatB <- plyr::ldply(diffmatB, rbind)
        diffmatB <- as.data.frame(diffmatB)
        
        if (dim(diffmatB)[1] > 0) {
            cnames <- colnames(diffmatB)
            cnames[1] <- "ISGroup"
            
            colnames(diffmatB) <- cnames
        }
        
        if (dim(diffmatB)[1] < dim(dataA)[1]) {
            diffmatC <- levelA_res2[-which(levelA_res$mz %in% diffmatB$mz), ]
            
            if (nrow(diffmatC) > 0) {
                isop_last <- paste0("ISgroup_", diffmatC$Module_RTclust, "_", 1)
                
                diffmatC <- cbind(isop_last, diffmatC)
                colnames(diffmatC) <- colnames(diffmatB)
                diffmatD <- rbind(diffmatB[, 1:10], diffmatC[, 1:10])
                
                rm(diffmatC)
                rm(diffmatB)
            } else {
                diffmatD <- diffmatB  #[,c(1:10)]
                rm(diffmatB)
            }
        } else {
            diffmatD <- diffmatB  #[,c(1:10)]
            rm(diffmatB)
        }
        rm(levelA_res2)
        diffmatD <- as.data.frame(diffmatD)
        diffmatD$mz <- as.numeric(diffmatD$mz)
        diffmatD <- diffmatD[order(diffmatD$mz, diffmatD$time), ]
        
        levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ]
        levelA_res1 <- cbind(diffmatD[, 1], levelA_res)
        
        isop_res_md <- cbind(diffmatD[, c(4, 5, 1, 3)], mean_int_vec, diffmatD[, 2])
        colnames(isop_res_md) <- c("mz", "time", "ISgroup", "Module_RTclust", "AvgIntensity", "MD")
        
        MD <- diffmatD[, 2]
        levelA_res1 <- cbind(levelA_res1[, 1:4], mean_int_vec, MD)
        rm(MD)
        
        cnames <- colnames(levelA_res1)
        cnames[1] <- "ISgroup"
        colnames(levelA_res1) <- cnames
        
        multiresmat <- merge(levelB_res, levelA_res1, by = "mz")
        multiresmat <- multiresmat[, c("mz", "time.x", "MatchCategory", "theoretical.mz", "chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "Module_RTclust", "time.y", "mean_int_vec", "MD")]
        colnames(multiresmat) <- c("mz", "time", "MatchCategory", "theoretical.mz", "chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "Module_RTclust", "time.y", "mean_int_vec", "MD")
        
        rm(levelB_res)
        
        multiresmat <- multiresmat[order(multiresmat$Module_RTclust), ]
        
        rm(m1)
        
        t2 <- table(multiresmat$mz)
        same1 <- which(t2 == 1)
        uniquemz <- names(same1)
        
        multiresmat$MatchCategory = rep("Multiple", dim(multiresmat)[1])
        multiresmat$MatchCategory[which(multiresmat$mz %in% uniquemz)] <- "Unique"
        
        setwd(outloc)
        
        multiresmat <- multiresmat[order(multiresmat$Module_RTclust, multiresmat$chemical_ID), ]
        dupmz <- multiresmat$mz[which(duplicated(multiresmat$mz) == TRUE)]
        tablemz <- table(multiresmat$Module_RTclust, multiresmat$mz)
        
        multiresmat$MatchCategory <- rep("Multiple", dim(multiresmat)[1])
        multiresmat$MatchCategory[-which(multiresmat$mz %in% dupmz)] <- "Unique"
        
        if (length(which(multiresmat$chemical_ID == "-")) > 0) {
            multiresmat <- multiresmat[-which(multiresmat$chemical_ID == "-"), ]
        }
        
        tall <- table(multiresmat$chemical_ID, multiresmat$mz)
        tall_checkmultimz <- apply(tall, 1, sum)
        tall_unimzperchem <- tall[which(tall_checkmultimz == 1), ]
        
        cnames <- colnames(dataA)
        cnames[2] <- "time"
        colnames(dataA) <- as.character(cnames)
        dataA <- unique(dataA)
        
        cnames <- colnames(multiresmat)
        cnames[10] <- "ISgroup"
        colnames(multiresmat) <- cnames
        
        level_module_isop_annot <- levelA_res1
        
        cnames <- colnames(multiresmat)
        cnames[2] <- "time"
        colnames(multiresmat) <- cnames
        
        rm(levelA_res)
        rm(levelB_res)
        rm(m2)
        
        mchemdata <- multiresmat
        
        rm(multiresmat)
        
        mchemdata <- as.data.frame(mchemdata)
        
        mchemdata$mz <- as.numeric(as.character(mchemdata$mz))
        mchemdata$time <- as.numeric(as.character(mchemdata$time))
        
        bad_rows <- which(abs(mchemdata$time - mchemdata$time.y) > 0)
        
        if (length(bad_rows) > 0) {
            mchemdata <- mchemdata[-bad_rows, ]
        }
        
        chemids <- unique(mchemdata$chemical_ID)
        chemids <- unique(chemids)
        
        if (length(chemids) > 10) {
            num_sets <- length(chemids)/2
        } else {
            num_sets <- 1
        }
        
        list_winsize <- num_sets
        list_size <- round(length(chemids)/list_winsize)
        
        if (length(chemids) > list_winsize) {
            g <- seq(1, length(chemids), list_size)
            g <- factor(g)
            chemids_split <- split(seq_along(chemids), f = g)
        } else {
            chemids_split <- split(seq_along(chemids), f = length(chemids))
        }
        
        num_sets <- length(chemids_split)
        
        mchemdata$mz <- round(as.numeric(as.character(mchemdata$mz)), 5)
        mchemdata$time <- round(as.numeric(as.character(mchemdata$time)), 1)
        mchemdata$MD <- round(as.numeric(as.character(mchemdata$MD)), 3)
        mchemdata$mean_int_vec <- round(as.numeric(as.character(mchemdata$mean_int_vec)), 1)
        mchemdata$time.y <- round(as.numeric(as.character(mchemdata$time.y)), 1)
        
        if (max_diff_rt >= 9999) {
            module_num <- gsub(mchemdata$Module_RTclust, pattern = "_[0-9]{1,}", replacement = "")
            module_num <- paste0(module_num, "_0")
            mchemdata$Module_RTclust <- module_num
        }
        write.csv(mchemdata, file = "Stage2.csv", row.names = FALSE)
        
        save(list = c("outloc", "adduct_weights", "boostIDs", "pathwaycheckmode", "adduct_table", "max_diff_rt", "corthresh", "filter.by", "max.rt.diff", "max_isp", "min_ions_perchem", "max.mz.diff", "db_name", "allsteps", "redundancy_check", "num_sets"), file = "tempobjects.Rda")
        save(list = c("mchemdata", "chemids", "adduct_table", "mzid", "max_diff_rt", "isop_res_md", "corthresh", "level_module_isop_annot", "chemids_split", "corthresh", "max.mz.diff", "outloc", "num_sets", "db_name", "num_nodes", "num_sets", "adduct_weights", "filter.by", 
            "max.rt.diff", "max_isp", "MplusH.abundance.ratio.check", "mass_defect_window", "mass_defect_mode", "allsteps"), file = "step1_results.Rda")
        
        rm(mchemdata)
        rm(chemids)
        rm(mzid)
        try(rm(global_cor), silent = TRUE)
        rm(isop_res_md)
        rm(level_module_isop_annot)
        rm(dataA)
        rm(tall_unimzperchem)
        rm(tablemz)
        rm(levelB_res2)
        rm(levelA_res1)
        rm(list = ls())
        load("tempobjects.Rda")
    } else {
        print("Status 1: Skipping step 1.")
        print("Status 2: Using existing step1_results.Rda file.")
        
        allsteps_temp <- allsteps
        load("tempobjects.Rda")
        allsteps <- allsteps_temp
    }
    
    if (allsteps == TRUE) {
        print("Status 3: Calculating scores for individual chemicals/metabolites")
        
        if (num_sets > num_nodes) {
            cl <- parallel::makeCluster(num_nodes)
            parallel::clusterEvalQ(cl, "multilevelannotationstep2")
            parallel::clusterExport(cl, "multilevelannotationstep2")
            parallel::clusterEvalQ(cl, "library(Rdisop)")
            parallel::clusterEvalQ(cl, "library(plyr)")
            parallel::clusterExport(cl, "get_chemscorev1.6.71")
            parallel::clusterExport(cl, "getMolecule")
            parallel::clusterExport(cl, "ldply")
            parallel::clusterExport(cl, "get_confidence_stage2")
            
            parallel::clusterExport(cl, "check_element")
            parallel::clusterExport(cl, "group_by_rt_histv2")
            parallel::clusterExport(cl, "adduct_table")
            parallel::clusterExport(cl, "adduct_weights")
            
            parallel::parLapply(cl, 1:num_sets, function(arg1) {
                cur_fname <- paste0(outloc, "/stage2/chem_score", arg1, ".Rda")
                check_if_exists <- suppressWarnings(try(load(cur_fname)))
                
                if (is(check_if_exists, "try-error")) {
                  xMSannotator::multilevelannotationstep2(outloc1 = outloc, list_number = arg1)
                } else {
                  print(paste0("List ", arg1, " already exists."))
                }
            })
            
            parallel::stopCluster(cl)
        } else {
            for (arg1 in 1:num_sets) {
                xMSannotator::multilevelannotationstep2(outloc1 = outloc, list_number = arg1)
            }
        }
        
        setwd(outloc)
        
        rm(list = ls())
        try(rm(hmdbCompMZ), silent = TRUE)
        
        load("tempobjects.Rda")
        
        print("Status 4: Pathway evaluation")
        multilevelannotationstep3(outloc = outloc, adduct_weights = adduct_weights, boostIDs = boostIDs, pathwaycheckmode = pathwaycheckmode)
        
        setwd(outloc)
        rm(list = ls())
        
        try(rm(hmdbCompMZ), silent = TRUE)
        try(rm(hmdbCompMZ, env = .GlobalEnv), silent = TRUE)
        
        load("tempobjects.Rda")
        
        print("Status 5: Assigning confidence levels")
        multilevelannotationstep4(outloc = outloc, max.rt.diff = max_diff_rt, filter.by = filter.by, adduct_weights = adduct_weights, min_ions_perchem = min_ions_perchem, boostIDs = boostIDs, max_isp = max_isp, max.mz.diff = max.mz.diff)
        
        rm(list = ls())
        
        load("tempobjects.Rda")
        
        try(rm(hmdbAllinf, env = .GlobalEnv), silent = TRUE)
        
        if (redundancy_check == TRUE) {
            print("Status 6:Redundancy filtering")
            
            rm(list = ls())
            
            load("tempobjects.Rda")
            
            suppressWarnings(multilevelannotationstep5(outloc = outloc, max.rt.diff = max_diff_rt, filter.by = filter.by, adduct_weights = adduct_weights, min_ions_perchem = min_ions_perchem, boostIDs = boostIDs, max_isp = max_isp, db_name = db_name, max.mz.diff = max.mz.diff))
        }
    }
    
    print("################")
    print("Final: Processing complete")
    print("Output files description:")
    print("Stage 1 includes clustering of features based on intensity and retention time without annotations")
    print("Stage 2 includes clustering results along with simple m/z based database matching")
    
    if (allsteps == TRUE) {
        print("Stage 3 includes scores for annotations assigned in stage 2 based on multiple criteria")
        print("Stages 4 and 5 include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
    }
    
    suppressWarnings(sink(file = NULL))
    
    setwd(outloc)
    sink(file = "Readme.txt")
    print("Output files description:")
    print("Stage1.csv includes clustering of features based on intensity and retention time without annotations")
    print("DBresults/Stage2.csv includes clustering results along with simple m/z based database matching")
    
    if (allsteps == TRUE) {
        print("DBresults/Stage3.csv includes scores for annotations assigned in stage 2 based on multiple criteria")
        print("DBresults/Stage4.csv and DBresults/Stage5.csv include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
    }
    suppressWarnings(sink(file = NULL))
}
