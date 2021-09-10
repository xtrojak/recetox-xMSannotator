get_confidence_stage4 <-function(curdata,
                                 max_diff_rt,
                                 adduct_weights = NA,
                                 filter.by = NA,
                                 min_ions_perchem = 1,
                                 max_isp = 5) {
    curdata <- curdata[order(curdata$Adduct), ]
    
    if (curdata$score[1] < 0.1) {
      score_level <- 0
      final_res <- cbind(score_level, curdata)
      final_res <- as.data.frame(final_res)
      return(final_res)
    }
    
    curdata_all <- curdata
    cur_adducts_with_isotopes <- curdata$Adduct
    
    cur_adducts <-gsub(cur_adducts_with_isotopes,
                       pattern = "(_\\[(\\+|\\-)[0-9]*\\])",
                       replacement = "")
    
    if (is.na(adduct_weights) == TRUE) {
      adduct_weights <- data.frame(Adduct = c("M+H", "M-H"), Weight = c(1, 1))
    }
  
    adduct_monoisot <- data.frame(cbind("M", 1, 1, 0, "-", "S"))
    colnames(adduct_monoisot) <- colnames(adduct_table)
    adduct_table <- rbind(adduct_table, adduct_monoisot)
    adduct_table <- adduct_table[order(adduct_table$Adduct), ]
    
    curdata <- cbind(curdata, cur_adducts)
    
    curdata <- merge(curdata, adduct_table, by.x = "cur_adducts", by.y = "Adduct")
    curdata <- curdata[, -c(1)]
    
    formula_vec <- curdata$Formula
    
    temp1 <- curdata
    temp1$Module_RTclust <- gsub(temp1$Module_RTclust, pattern = "(_[0-9]*)", replacement = "")
    
    table_modules <- table(temp1$Module_RTclust)
    
    rm(temp1)
    module_names <- names(table_modules[which(table_modules > 0)])
    chemscoremat_conf_levels <- "High"
    if (length(which(adduct_table$Adduct %in% cur_adducts)) < 1)
    {
      chemscoremat_conf_levels <- "None"
    }
    if (length(unique(curdata$Adduct)) < 2 &&
        curdata$score < 10 && length(which(cur_adducts %in% filter.by)) < 1) {
      chemscoremat_conf_levels <- "None"
    } else{
      if (length(unique(curdata$Adduct)) < 2 &&
          length(which(cur_adducts %in% filter.by)) > 0) {
        chemscoremat_conf_levels <- "Low"
      }
    }
    
    if (length(unique(curdata$Adduct)) < 2 && nrow(curdata) > 1) {
      if (nrow(curdata) == 2) {
        chemscoremat_conf_levels <- "Medium"
      }
      
      if (length(which(cur_adducts %in% filter.by)) < 1)
      {
        chemscoremat_conf_levels <- "None"
      } else {
        curdata <- curdata[which(cur_adducts %in% filter.by), ]
        score_level <- "2"
        score_level <- as.numeric(as.character(score_level))
        final_res <- cbind(score_level, curdata)
        return(final_res)
      }
    }
    
    curdata$time <- as.numeric(as.character(curdata$time))
    delta_rt <- max(curdata$time) - min(curdata$time)
    
    delta_rt <- round(delta_rt)
    
    if (delta_rt > max_diff_rt)
    {
      chemscoremat_conf_levels <- "None"
      
      groupB <-group_by_rt(curdata, time_step = 3, max.rt.diff = max_diff_rt, groupnum = "")
      groupB <- groupB[order(groupB$mz), ]
      
      cname1 <- paste(groupB$mz, groupB$time, sep = "_")
      if (length(which(duplicated(cname1) == TRUE)) > 0) {
        groupB <- groupB[-which(duplicated(cname1) == TRUE), ]
      }
      curdata <- curdata[order(curdata$mz), ]

      module_clust <- groupB[, 1]
      curdata$Module_RTclust <- module_clust
      
      if (length(which(adduct_weights[, 1] %in% cur_adducts)) > 0 && curdata$score[1] > 0.1) {
        if (is.na(filter.by) == TRUE) {
          good_mod <- curdata$Module_RTclust[which(curdata$Adduct %in% adduct_weights[, 1])]
          curdata <- curdata[which(curdata$Module_RTclust %in% good_mod), ]
          t1 <- table(curdata$Module_RTclust)
          t2 <- t1[which(t1 > 0)]
          t2 <- t2[which(t2 == max(t2)[1])]
          n1 <- names(t2)
          curdata <- curdata[which(curdata$Module_RTclust == n1[1]), ]
          delta_rt <- max(curdata$time) - min(curdata$time)
          delta_rt <- round(delta_rt)
          if (curdata$score[1] > 0 &&
              nrow(curdata) > 1 &&
              length(unique(curdata$Adduct)) > 1 && delta_rt < max_diff_rt) {
            chemscoremat_conf_levels <- "High"
          } else {
            chemscoremat_conf_levels <- "Low"
          }
        } else {
          if (length(which(cur_adducts %in% filter.by)) > 0) {
            good_mod <- curdata$Module_RTclust[which(curdata$Adduct %in% filter.by)]
            curdata <- curdata[which(curdata$Module_RTclust %in% good_mod), ]
            curdata <- as.data.frame(curdata)
            t1 <- table(curdata$Module_RTclust)
            t2 <- t1[which(t1 > 0)]
            t2 <- t2[which(t2 == max(t2)[1])]
            n1 <- names(t2)
            curdata <- curdata[which(curdata$Module_RTclust == n1[1]), ]
            delta_rt <- max(curdata$time) - min(curdata$time)
            delta_rt <- round(delta_rt)
            if (length(which(t1 > 1)) < 1) {
              if (curdata$score[1] > 0 &&
                  nrow(curdata) > 1 &&
                  length(unique(curdata$Adduct)) > 1 &&
                  delta_rt < max_diff_rt &&
                  length(which(curdata$Adduct %in% filter.by)) > 0) {
                chemscoremat_conf_levels <- "High"
              } else {
                if (curdata$score[1] > 0 &&
                    length(which(curdata$Adduct %in% filter.by)) > 0) {
                  chemscoremat_conf_levels <- "Medium"
                } else {
                  if (curdata$score[1] > 0) {
                    chemscoremat_conf_levels <- "Low"
                    
                  } else {
                    chemscoremat_conf_levels <- "None"
                  }
                }
              }
            } else {
              t2 <- t1[which(t1 > 0)]
              t2 <- t2[which(t2 == max(t2)[1])]
              n1 <- names(t2)
              curdata <- curdata[which(curdata$Module_RTclust == n1[1]), ]
              delta_rt <- max(curdata$time) - min(curdata$time)
              delta_rt <- round(delta_rt)
              if (curdata$score[1] > 0 &&
                  nrow(curdata) > 1 &&
                  length(unique(curdata$Adduct)) > 1 &&
                  delta_rt < max_diff_rt &&
                  length(which(curdata$Adduct %in% filter.by)) > 0) {
                chemscoremat_conf_levels <- "High"
              } else {
                if (curdata$score[1] > 0 &&
                    length(which(curdata$Adduct %in% filter.by)) > 0 &&
                    delta_rt < max_diff_rt) {
                  chemscoremat_conf_levels <- "Medium"
                } else {
                  if (curdata$score[1] > 0) {
                    chemscoremat_conf_levels <- "Low"
                  } else {
                    chemscoremat_conf_levels <- "None"
                  }
                }
              }
            }
          }
        }
      }
      else {
        chemscoremat_conf_levels <- "None"
      }
      module_clust <- paste(curdata$module_num, curdata$Module_RTclust, sep = "")
      curdata$Module_RTclust <- module_clust
    }
    curdata <- as.data.frame(curdata)
    
    temp1 <- curdata
    temp1$Module_RTclust <- gsub(temp1$Module_RTclust, pattern = "(_[0-9]*)", replacement = "")
    table_modules <- table(temp1$Module_RTclust)
    
    rm(temp1)
    module_names <- names(table_modules[which(table_modules > 0)])
    {
      #change in 1.5.3
      if (length(module_names) > 1 && curdata$score < 10) {
        chemscoremat_conf_levels <- "None"
      } else {
        if (length(module_names) > 1 &&
            curdata$score > 10 && length(which(cur_adducts %in% filter.by)) > 0) {
          chemscoremat_conf_levels <- "Medium"
        }
      }
      {
        if (nrow(curdata) < 2) {
          if (length(which(cur_adducts %in% adduct_weights[, 1])) < 1) {
            chemscoremat_conf_levels <- "Low"
          } else {
            # matches an M+H and score is greater than 10
            if (curdata$score > 10 &&
                length(which(cur_adducts %in% filter.by)) > 0) {
              chemscoremat_conf_levels <- "Medium"
            } else {
              chemscoremat_conf_levels <- "None"
            }
          }
        } else {
          min_nummol <- min(curdata$num_molecules)
          min_nummol_ind <- which(curdata$num_molecules == min_nummol)
          # num molecules check
          if (min_nummol > 1) {
            chemscoremat_conf_levels <- "Low"
          } else {
            # Multiple molecules abundance check
            check1 <- gregexpr(text = cur_adducts, pattern = "([2-3]+M)")
            if (length(check1) > 0) {
              min_mol <- min(curdata$num_molecules)
              min_mol_ind <- which(curdata$num_molecules == min_mol)
              
              max_int_min_mol <- max(curdata$mean_int_vec[min_mol_ind])
              
              bad_ind_status <- rep(0, length(check1))
              
              min_mol <- min(adduct_table[which(adduct_table$Adduct %in% cur_adducts), 2])
              
              if (min_mol < 2) {
                # check pattern of each adduct
                for (a1 in 1:length(check1)) {
                  strlength <- attr(check1[[a1]], "match.length")
                  if (strlength[1] > (-1)) {
                    abund_ratio <- curdata$mean_int_vec[a1] / max_int_min_mol
                    
                    if (is.na(abund_ratio) == FALSE) {
                      if (abund_ratio > 1) {
                        bad_ind_status[a1] <- 1
                        if (chemscoremat_conf_levels == "High") {
                          chemscoremat_conf_levels <- "Medium"
                        } else {
                          chemscoremat_conf_levels <- "Low"
                        }
                      }
                    } else {
                      chemscoremat_conf_levels <- "Low"
                      bad_ind_status[a1] <- 1
                    }
                  }
                }
                
                if (length(which(bad_ind_status == 1)) > 0) {
                  bad_adducts <- cur_adducts[which(bad_ind_status == 1)]
                }
                
                if (length(nrow(curdata) > 0)) {
                  cur_adducts_with_isotopes <- curdata$Adduct
                  
                  cur_adducts <-gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
                } else {
                  chemscoremat_conf_levels <- "None"
                }
              } else {
                chemscoremat_conf_levels <- "Low"
              }
              adduct_charges <- curdata$charge
              min_charge <- min(curdata$charge)
              min_charge_ind <- which(curdata$charge == min_charge)
              
              high_charge_ind <- which(curdata$charge > 1)
              max_int_min_charge <- max(curdata$mean_int_vec)
              
              if (length(high_charge_ind) > 0) {
                abund_ratio_min_maxcharge <-
                  max(curdata$mean_int_vec[min_charge_ind])[1] / max(curdata$mean_int_vec[high_charge_ind])[1]
                
                if ((abund_ratio_min_maxcharge < 1))
                {
                  chemscoremat_conf_levels <- "Low"
                }
              }
              
              bad_ind_status <- rep(1, length(check1))
              
              if ((min_charge > 1))
              {
                chemscoremat_conf_levels <- "Low"
              }
              
              # isotope based check for +1 and +2 to make sure that the abundant form is present
              check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
              
              if (length(check2) > 0) {
                for (a1 in 1:length(check2)) {
                  strlength <- attr(check2[[a1]], "match.length")
                  
                  if (strlength[1] > (-1)) {
                    count_abundant_form <- length(which(cur_adducts %in% cur_adducts[a1]))
                    
                    if (count_abundant_form < 2)
                    {
                      curdata <- curdata[-c(a1), ]
                    }
                  }
                }
              }
              
              cur_adducts_with_isotopes <- curdata$Adduct

              cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
              formula_vec <- curdata$Formula
              curformula <- as.character(formula_vec[1])
            }
          }
        }
      }
    }
    formula_vec <- curdata$Formula
    curformula <- as.character(formula_vec[1])
    curformula <- gsub(curformula, pattern = "Ca", replacement = "")
    curformula <- gsub(curformula, pattern = "Cl", replacement = "")
    curformula <- gsub(curformula, pattern = "Cd", replacement = "")
    curformula <- gsub(curformula, pattern = "Cr", replacement = "")
    
    numcarbons <- check_element(curformula, "C")
    if (numcarbons < 1) {
      chemscoremat_conf_levels <- "None"
    }
    
    if (length(unique(curdata$Adduct)) < 2 && nrow(curdata) > 1) {
      if (nrow(curdata) == 2) {
        chemscoremat_conf_levels <- "Medium"
      }
      
      if (length(which(cur_adducts %in% filter.by)) < 1)
      {
        chemscoremat_conf_levels <- "None"
        return(chemscoremat_conf_levels)
      }
    } else {
      if (length(unique(curdata$Adduct)) < 2 &&
          curdata$score < 10 && length(which(cur_adducts %in% filter.by)) < 1) {
        chemscoremat_conf_levels <- "None"
      } else {
        if (length(unique(curdata$Adduct)) < 2 &&
            curdata$score > 10 && length(which(cur_adducts %in% filter.by)) > 0) {
          chemscoremat_conf_levels <- "Medium"
        } else {
          if (length(unique(curdata$Adduct)) < 2 &&
              curdata$score < 10 && length(which(cur_adducts %in% filter.by)) > 0) {
            chemscoremat_conf_levels <- "Low"
          }
        }
      }
    }
    
    if (nrow(curdata) < 1) {
      score_level <- 0
      curdata <- curdata_all
    } else {
      # 3: High
      # 2: medium
      # 1: Low
      # 0: None
      
      score_level <- as.character(chemscoremat_conf_levels)
      
      score_level <- gsub(score_level, pattern = "High", replacement = "3")
      score_level <- gsub(score_level, pattern = "Medium", replacement = "2")
      score_level <- gsub(score_level, pattern = "Low", replacement = "1")
      score_level <- gsub(score_level, pattern = "None", replacement = "0")
    }
    
    if (is.na(score_level[1]) == TRUE) {
      stop(curdata)
    }
    
    score_level <- as.numeric(as.character(score_level))
    final_res <- cbind(score_level, curdata)
    
    final_res <- as.data.frame(final_res)
    
    if (length(which(is.na(final_res$score_level) == TRUE)) > 0) {
      final_res$score_level <-
        replace(final_res$score_level, which(is.na(final_res$score_level) == TRUE), 0)
    }
    
    num_uniq_adducts <- length(unique(curdata$Adduct))
    uniq_adducts <- unique(curdata$Adduct)
    num_good_adducts <- length(which(uniq_adducts %in% filter.by))
    
    if (num_good_adducts > 0) {
      final_res$score <- final_res$score * num_good_adducts
    } else {
      final_res$score <- final_res$score * (0)
    }
    
    if (nrow(final_res) < min_ions_perchem) {
      final_res$score_level <- 0
    }

    final_res <- as.data.frame(final_res)
    
    rm(
      list = c(
        "curdata_all",
        "cur_adducts_with_isotopes",
        "cur_adducts",
        "adduct_monoisot",
        "curdata",
        "temp1",
        "table_modules",
        "chemscoremat_conf_levels",
        "score_level",
        "groupB",
        "good_mod",
        "module_clust",
        "uniq_adducts",
        "num_good_adducts",
        "num_uniq_adducts"
      )
    )
    
    return(final_res)
  }
