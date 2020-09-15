get_chemscorev1.6.71 <- function(chemicalid, mchemicaldata, corthresh, global_cor, mzid, max_diff_rt, level_module_isop_annot, adduct_weights, filter.by, max_isp, MplusH.abundance.ratio.check, mass_defect_window, mass_defect_mode) {
  mchemicaldata$mz <- as.numeric(as.character(mchemicaldata$mz))
  mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))
  mchemicaldata$MonoisotopicMass <- as.numeric(as.character(mchemicaldata$MonoisotopicMass))

  level_module_isop_annot$mz <- as.numeric(as.character(level_module_isop_annot$mz))
  level_module_isop_annot$time <- as.numeric(as.character(level_module_isop_annot$time))

  mchemicaldata$Module_RTclust <- gsub(mchemicaldata$Module_RTclust, pattern = "_[0-9]*", replacement = "")

  if (length(which(mchemicaldata$Adduct %in% as.character(adduct_weights[, 1]))) > 0) {
    rel_adduct_data <- mchemicaldata[which(mchemicaldata$Adduct %in% as.character(adduct_weights[, 1])),]
    rel_adduct_module <- gsub(rel_adduct_data$Module_RTclust, pattern = "_[0-9]*", replacement = "")
    module_rt_group <- gsub(mchemicaldata$Module_RTclust, pattern = "_[0-9]*", replacement = "")

    mchemicaldata <- mchemicaldata[which(module_rt_group %in% rel_adduct_module),]
  }

  mchemicaldata <- unique(mchemicaldata)

  curformula <- as.character(mchemicaldata[, 7])
  formula_check <- Rdisop::getMolecule(as.character(curformula[1]))
  exp_isp <- which(formula_check$isotopes[[1]][2,] >= 0.001)
  abund_ratio_vec <- formula_check$isotopes[[1]][2, exp_isp]

  numoxygen <- check_element(curformula, "O")
  water_adducts <- c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
  water_adduct_ind <- which(mchemicaldata$Adduct %in% water_adducts)

  if (numoxygen < 1) {
    if (length(water_adduct_ind) > 0) {
      mchemicaldata <- mchemicaldata[-water_adduct_ind,]
    }
  }

  #water check
  if (nrow(mchemicaldata) < 1) {
    return(list("chemical_score" = -100, "filtdata" = mchemicaldata))
  }
  mchemicaldata$Adduct <- as.character(mchemicaldata$Adduct)

  final_isp_annot_res_all <- mchemicaldata
  level_module_isop_annot <- as.data.frame(level_module_isop_annot)
  level_module_isop_annot$Module_RTclust <- gsub(level_module_isop_annot$Module_RTclust, pattern = "_[0-9]*", replacement = "")

  mchemicaldata_goodadducts_index <- which(mchemicaldata$Adduct %in% as.character(adduct_weights[, 1]))

  if (length(mchemicaldata_goodadducts_index) > 0) {
    final_isp_annot_res_isp <- lapply(seq_along(mchemicaldata_goodadducts_index), function(i) {
      m <- mchemicaldata_goodadducts_index[i]
      final_isp_annot_res <- cbind(paste0("group", i), mchemicaldata[m,])
      module_rt_group <- as.character(mchemicaldata$Module_RTclust[m])
      module_rt_group <- gsub(module_rt_group, pattern = "_[0-9]*", replacement = "")

      isp_mat_module_rt_group <- as.character(level_module_isop_annot$Module_RTclust)

      query_md <- mchemicaldata$mz[m] - round(mchemicaldata$mz[m])
      query_rt <- mchemicaldata$time[m]
      query_int <- mchemicaldata$mean_int_vec[m]

      put_isp_masses_curmz_data <- level_module_isop_annot[
        which(abs(level_module_isop_annot$time - query_rt) < max_diff_rt &
          abs((level_module_isop_annot$MD) - (query_md)) < mass_defect_window &
          isp_mat_module_rt_group == module_rt_group &
          level_module_isop_annot$AvgIntensity < query_int),]

      put_isp_masses_curmz_data <- as.data.frame(put_isp_masses_curmz_data)
      put_isp_masses_curmz_data$mz <- as.numeric(as.character(put_isp_masses_curmz_data$mz))
      put_isp_masses_curmz_data$time <- as.numeric(as.character(put_isp_masses_curmz_data$time))
      put_isp_masses_curmz_data <- unique(put_isp_masses_curmz_data)
      put_isp_masses_curmz_data <- put_isp_masses_curmz_data[order(put_isp_masses_curmz_data$mz),]

      if (length(put_isp_masses_curmz_data) > 0) {
        int_vec <- put_isp_masses_curmz_data[, 5] / mchemicaldata[m, 13]
        put_isp_masses_curmz_data[, 2] <- as.numeric(as.character(put_isp_masses_curmz_data[, 2]))

        ischeck <- which(int_vec <= max(abund_ratio_vec[-1] + 0.10))

        for (isp_v in ischeck) {
          diff_rt <- abs(put_isp_masses_curmz_data[isp_v, 2] - mchemicaldata[m, 2])
          isnum <- (round(put_isp_masses_curmz_data[isp_v, 1]) - round(mchemicaldata[m, 1]))

          if (mass_defect_mode == "neg" | mass_defect_mode == "both") {
            isnum <- abs(isnum)
          }

          if (diff_rt < max_diff_rt &
            isnum <= max_isp &
            max_isp > 0 &
            isnum <= max_isp &
            isnum > 0) {
            isnum2 <- round(put_isp_masses_curmz_data[isp_v, 1]) - round(mchemicaldata$MonoisotopicMass[m])
            isnum2 <- round(isnum2)

            isp_sign <- if (isnum2 > 0) "+" else "-"

            isnum2 <- abs(isnum2)

            if (isnum2 <= max_isp) {
              form_name <- as.character(paste0(mchemicaldata[m, 7], "_[", isp_sign, (isnum2), "]"))

              other_inf <- cbind(rep("-", 7))
              temp_var <- cbind(put_isp_masses_curmz_data[isp_v, 1:2], t(other_inf), put_isp_masses_curmz_data[isp_v, 3:4], put_isp_masses_curmz_data[isp_v, c(2, 5:6)])
              temp_var <- as.data.frame(temp_var)
              colnames(temp_var) <- colnames(mchemicaldata)

              temp_var$Formula <- form_name
              temp_var$Name <- as.character(mchemicaldata[m, 6])
              temp_var$chemical_ID <- as.character(mchemicaldata[m, 5])
              temp_var$Adduct <- paste0(mchemicaldata[m, 9], "_[", isp_sign, (abs(isnum)), "]")
              temp_var <- as.data.frame(temp_var)
              temp_var <- cbind(paste0("group", i), temp_var)

              final_isp_annot_res <- as.data.frame(final_isp_annot_res)

              if (nrow(temp_var) > 0) {
                check_mz <- which(temp_var$mz %in% final_isp_annot_res)

                if (length(check_mz) > 0) {
                  temp_var <- temp_var[-check_mz,]
                }
                if (nrow(temp_var) > 0) {
                  final_isp_annot_res <- rbind(final_isp_annot_res, temp_var)
                }
              }
            }
          }
        }
      }
      return(final_isp_annot_res)
    })

    rm(level_module_isop_annot)

    final_isp_annot_res2 <- ldply(final_isp_annot_res_isp, rbind)
    final_isp_annot_res2 <- as.data.frame(final_isp_annot_res2)
    final_isp_annot_res2 <- final_isp_annot_res2[, -1]
    final_isp_annot_res2 <- as.data.frame(final_isp_annot_res2)

    mchemicaldata <- rbind(final_isp_annot_res_all, final_isp_annot_res2)
  }
  mchemicaldata <- unique(mchemicaldata)

  bad_rows <- which(is.na(mchemicaldata$mz) == TRUE)
  if (length(bad_rows) > 0) {
    mchemicaldata <- mchemicaldata[-bad_rows,]
  }
  mchemicaldata <- mchemicaldata[order(mchemicaldata$mz),]

  mod_names <- mchemicaldata$Module_RTclust
  mod_names <- unique(mod_names)

  diffmatB <- lapply(mod_names, function(groupA_num) {
    subdata <- mchemicaldata[which(mchemicaldata$Module_RTclust == groupA_num),]
    subdata <- subdata[order(subdata$time),]

    if (nrow(subdata) > 0) {
      subdata <- group_by_rt_histv2(subdata, time_step = 1, max_diff_rt = 10, groupnum = groupA_num)
    }
    rownames(subdata) <- NULL

    return(subdata)
  })

  mchemicaldata <- ldply(diffmatB, rbind)
  mchemicaldata <- unique(mchemicaldata)

  rm(diffmatB)
  rm(final_isp_annot_res)

  #write.table(mchemicaldata, file = "../Stage2_withisotopes.txt", append = TRUE, sep = "\t", col.names = FALSE)

  table_mod <- table(mchemicaldata$Module_RTclust)
  table_mod <- table_mod[table_mod > 0]
  table_mod <- table_mod[order(table_mod, decreasing = TRUE)]

  mchemicaldata_orig <- mchemicaldata
  top_mod <- names(table_mod)

  best_chemical_score <- -100

  for (i in seq_along(which(table_mod >= 1))) {
    dup_add <- { }
    chemical_score <- -99999
    mchemicaldata <- mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust == top_mod[i]),]
    mchemicaldata <- mchemicaldata[order(mchemicaldata$mz),]

    cur_adducts_with_isotopes <- mchemicaldata$Adduct

    if (nrow(mchemicaldata) < 2) {
      good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))

      if (good_adducts_len > 0) {
        chemical_score <- 10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))
        chemical_score <- chemical_score[1]

        if (chemical_score > best_chemical_score) {
          best_chemical_score <- chemical_score
          best_data <- mchemicaldata
        } else if (chemical_score == best_chemical_score) {
          best_chemical_score <- chemical_score
          best_data <- rbind(best_data, mchemicaldata)
        }
      }
    }

    check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
    mzid_cur <- paste(mchemicaldata$mz, mchemicaldata$time, sep = "_")

    dup_features <- which(duplicated(mzid_cur) == TRUE)
    if (length(dup_features) > 0) {
      mchemicaldata <- mchemicaldata[-dup_features,]
      mzid_cur <- paste(mchemicaldata$mz, mchemicaldata$time, sep = "_")
    }

    cor_mz <- global_cor[which(mzid %in% mzid_cur), which(mzid %in% mzid_cur)]
    cor_mz <- round(cor_mz, 1)

    if (length(cor_mz) > 1) {
      corrownamesA <- rownames(cor_mz)
      mat_rownames <- strsplit(as.character(corrownamesA), split = "_")

      m1 <- { }
      for (i in seq_along(mat_rownames)) {
        m1 <- rbind(m1, cbind(mat_rownames[[i]][1], mat_rownames[[i]][2]))
      }

      m1 <- as.data.frame(m1)
      colnames(m1) <- c("mz", "time")
      cor_mz2 <- cbind(m1, cor_mz)

      mz_order <- order(cor_mz2$mz)
      cor_mz2 <- cor_mz2[mz_order,]
      cor_mz <- cor_mz2[, -(1:2)]
      cor_mz <- cor_mz[, mz_order]
    }

    if (length(cor_mz) > 1) {
      check_cor <- sapply(1:dim(cor_mz)[1], function(k) {
        count_mz <- length(which(cor_mz[k,] >= corthresh)) - 1
        return(count_mz)
      })

      topquant_cor <- max(cor_mz[upper.tri(cor_mz)])
      check_cor2 <- check_cor / length(check_cor)

      if (length((cur_adducts_with_isotopes %in% adduct_weights[, 1] == TRUE)) > 0) {
        check_cor[check2 > 0 | (cur_adducts_with_isotopes %in% adduct_weights[, 1] == FALSE)] <- 0
      }

      #at least score of 2
      if (length(which(check_cor > 0) == TRUE) > 0) {
        hub_mz_list <- which(check_cor > 0 & (cur_adducts_with_isotopes %in% filter.by == TRUE))

        if (length(hub_mz_list) < 1) {
          hub_mz_list <- which(check_cor > 0 &
            check2 < 0 &
            (cur_adducts_with_isotopes %in% adduct_weights[, 1] == TRUE))
        }

        if (length(hub_mz_list) < 1) {
          hub_mz_list <- which(check_cor > 0 &
            check2 < 0 &
            (cur_adducts_with_isotopes %in% adduct_weights[, 1] == TRUE))
        }

        if (length(hub_mz_list) < 1) {
          hub_mz_list <- which(check_cor > 0 & check2 < 0)
        }

        hub_mz_int <- hub_mz_list[which(mchemicaldata$mean_int_vec[hub_mz_list] == max(mchemicaldata$mean_int_vec[hub_mz_list]))[1]]

        max_time_neighbors <- 0
        best_hub_time_mz <- hub_mz_int

        for (h1 in hub_mz_list) {
          hub_rt <- mchemicaldata$time[h1]

          diff_rt_hubmz <- apply(mchemicaldata, 1, function(k) {
            curtime <- as.numeric(as.character(k[3]))
            return(abs(hub_rt - curtime))
          })
          num_time_neighbors <- length(which(diff_rt_hubmz <= max_diff_rt))

          if (num_time_neighbors > max_time_neighbors) {
            best_hub_time_mz <- h1
            max_time_neighbors <- num_time_neighbors
          }
        }

        hub_mz <- best_hub_time_mz
        hub_rt <- mchemicaldata$time[hub_mz]

        diff_rt_hubmz <- apply(mchemicaldata, 1, function(k) {
          curtime <- as.numeric(as.character(k[3]))
          return(abs(hub_rt - curtime))
        })

        if (MplusH.abundance.ratio.check == TRUE) {
          layer_one_associations <- which(cor_mz[hub_mz,] >= corthresh &
            mchemicaldata$mean_int_vec < mchemicaldata$mean_int_vec[hub_mz] &
            diff_rt_hubmz <= max_diff_rt)
        } else {
          layer_one_associations <- which(cor_mz[hub_mz,] >= corthresh & diff_rt_hubmz <= max_diff_rt)
        }

        selected_mz <- c(hub_mz, layer_one_associations)
        selected_mz <- unique(selected_mz)

        if (is.na(topquant_cor) == TRUE) {
          topquant_cor <- 0
        }

        if (length(which(check_cor2 >= 0.1)) > 0) {
          mchemicaldata <- mchemicaldata[selected_mz,]
        }else {
          mchemicaldata <- mchemicaldata[which(cor_mz[hub_mz,] >= corthresh),]
        }

        mchemicaldata <- na.omit(mchemicaldata)
        if (nrow(mchemicaldata) < 2) {
          next
        }
        mchemicaldata <- as.data.frame(mchemicaldata)

        diff_rt <- abs(min(as.numeric(mchemicaldata$time)) - max(as.numeric(mchemicaldata$time)))
        diff_rt <- round(diff_rt)

        if (length(which(is.na(mchemicaldata$time)) == TRUE) > 0) {
          mchemicaldata <- mchemicaldata[-which(is.na(mchemicaldata$time) == TRUE),]
        }

        if (nrow(mchemicaldata) < 2) {
          next
        }

        if (diff_rt <= max_diff_rt) {
          dup_add <- which(duplicated(mchemicaldata$Adduct))
          if (length(dup_add) > 0) {
            dup_data <- mchemicaldata[dup_add,]
          }

          if (dim(mchemicaldata)[1] > 1) {
            mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))

            cur_adducts_with_isotopes <- mchemicaldata$Adduct
            cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])", replacement = "")

            good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))

            #change
            chemical_score <- length(unique(cur_adducts)) *
              good_adducts_len *
              (1 * (topquant_cor))

            if (good_adducts_len > 0) {
              chemical_score <- sum(chemical_score * (10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))))
              chemical_score <- chemical_score[1]
            }

            check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]\\])")

            if (length(which(check2 > 0)) > 0) {
              chemical_score <- 100 * chemical_score
            }
            names(chemical_score) <- chemicalid[1]
          }
        } else {
          mchemicaldata$Module_RTclust <- gsub(mchemicaldata$Module_RTclust, pattern = "_[0-9]*", replacement = "")
          mchemicaldata <- cbind(mchemicaldata[, 2:11], mchemicaldata[, 1], mchemicaldata[, 12:14])
          colnames(mchemicaldata) <- c("mz", "time", "MatchCategory", "theoretical.mz", "chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "Module_RTclust", "time.y", "mean_int_vec", "MD")

          mchemicaldata <- as.data.frame(mchemicaldata)
          mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))

          groupnumA <- unique(mchemicaldata$Module_RTclust)

          mchemicaldata <- group_by_rt_histv2(mchemicaldata, time_step = 1, max_diff_rt = max_diff_rt, groupnum = groupnumA)

          top_mod_sub <- table(mchemicaldata$Module_RTclust)
          top_mod_sub_names <- names(top_mod_sub)
          max_top_mod <- which.max(top_mod_sub)

          mchemicaldata <- mchemicaldata[which(mchemicaldata$Module_RTclust == top_mod_sub_names[max_top_mod]),]
          mchemicaldata <- mchemicaldata[order(mchemicaldata$mz),]

          s1 <- summary(mchemicaldata$time)
          iqr1 <- s1[5] - s1[2]
          min_val <- s1[2] - (1.5 * iqr1)
          max_val <- s1[5] + (1.5 * iqr1)

          if (min_val < s1[1] && max_val > s1[6]) {
            iqr1 <- min(abs(s1[3] - s1[2]), abs(s1[3] - s1[5]))
            min_val <- s1[2] - (1.5 * iqr1)
            max_val <- s1[5] + (1.5 * iqr1)
          }

          if (min_val < s1[1]) {
            min_val <- s1[1]
          }
          if (max_val > s1[6]) {
            max_val <- s1[6]
          }

          iqr1 <- min(abs(s1[3] - s1[2]), abs(s1[3] - s1[5]))
          iqr1 <- max(iqr1, max_diff_rt)

          diff_rt <- abs(max(mchemicaldata$time) - min(mchemicaldata$time))

          if (nrow(mchemicaldata) < 1) {
            next
          }

          if (diff_rt > 2 * max_diff_rt) {
            time_cor_groups <- sapply(list(myData1 = mchemicaldata), function(x)  split(x, cut(mchemicaldata$time, breaks = seq(min_val - max_diff_rt, max_val + max_diff_rt, iqr1))))
          } else {
            max_val <- max(mchemicaldata$time)
            min_val <- min(mchemicaldata$time)
            diff_rt <- abs(max(mchemicaldata$time) - min(mchemicaldata$time))

            if (min_val < max_diff_rt) {
              time_cor_groups <- sapply(list(myData1 = mchemicaldata), function(x) split(x, cut(mchemicaldata$time, breaks = c(0, max_val + 1))))
            } else if (diff_rt < max_diff_rt) {
              time_cor_groups <- sapply(list(myData1 = mchemicaldata), function(x) split(x, cut(mchemicaldata$time, breaks = c(min_val - diff_rt, max_val + diff_rt, diff_rt * 2))))
            } else {
              time_cor_groups <- sapply(list(myData1 = mchemicaldata), function(x)  split(x, cut(mchemicaldata$time, breaks = seq(min_val - diff_rt, max_val + diff_rt, 1 * max_diff_rt))))
            }
          }

          group_sizes <- sapply(time_cor_groups, function(x) nrow(as.data.frame(x)))

          if (max(group_sizes) < 2) {
            mchemicaldata <- mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust == top_mod[i]),]
            mchemicaldata <- mchemicaldata[order(mchemicaldata$mz),]

            cur_adducts_with_isotopes <- mchemicaldata$Adduct
            cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")

            good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))

            #change
            chemical_score <- length(unique(cur_adducts)) *
              good_adducts_len *
              (1 * (topquant_cor))

            if (good_adducts_len > 0) {
              chemical_score <- sum(chemical_score * (as.numeric(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts), 2])))
            }
          }

          if (any(group_sizes > 1)) {
            temp_best_score <- -100
            temp_best_data <- { }

            good_temp <- which(sapply(time_cor_groups, function(x) {
              any(x$Adduct %in% adduct_weights[, 1])
            }))

            for (mchemicaldata in time_cor_groups[good_temp]) {
              if (length(mchemicaldata) < 1) {
                next
              }

              if (nrow(mchemicaldata) < 1) {
                next
              }

              diff_rt <- abs(min(as.numeric(mchemicaldata$time)) - max(as.numeric(mchemicaldata$time)))
              diff_rt <- round(diff_rt)

              if (diff_rt <= max_diff_rt) {
                if (nrow(mchemicaldata) > 1) {
                  mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))
                  cur_adducts_with_isotopes <- mchemicaldata$Adduct
                  cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])", replacement = "")

                  good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
                  chemical_score <- length(unique(cur_adducts)) *
                    good_adducts_len *
                    (1 * (topquant_cor))

                  if (good_adducts_len > 0) {
                    chemical_score <- sum(chemical_score * (10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))))
                    chemical_score <- chemical_score[1]
                  }

                  check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])")

                  if (length(which(check2 > 0)) > 0) {
                    chemical_score <- 100 * chemical_score
                  }
                } else {
                  chemical_score <- 0
                }
              } else {
                d1 <- density(mchemicaldata$time, bw = max_diff_rt, from = min(mchemicaldata$time) - 0.001, to = (0.01 + max(mchemicaldata$time)), na.rm = TRUE)
                s1 <- summary(d1$x)

                min_val <- s1[2]
                max_val <- s1[5]

                if (length(which(mchemicaldata$time >= min_val & mchemicaldata$time <= max_val)) > 1) {
                  mchemicaldata <- mchemicaldata[which(mchemicaldata$time >= (min_val - 1) & mchemicaldata$time <= (max_val - 1)),]

                  if (nrow(mchemicaldata) > 1) {
                    mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))
                    cur_adducts_with_isotopes <- mchemicaldata$Adduct
                    cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])", replacement = "")

                    good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
                    chemical_score <- length(unique(cur_adducts)) *
                      good_adducts_len *
                      (1 * (topquant_cor))

                    if (good_adducts_len > 0) {
                      chemical_score <- sum(chemical_score * (10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))))
                      chemical_score <- chemical_score[1]
                    }

                    cur_adducts_with_isotopes <- mchemicaldata$Adduct
                    check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]\\])")

                    if (length(which(check2 > 0)) > 0) {
                      chemical_score <- 100 * chemical_score
                    }
                  }
                } else {
                  cur_adducts_with_isotopes <- mchemicaldata$Adduct
                  cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])", replacement = "")

                  good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
                  chemical_score <- length(unique(cur_adducts)) *
                    good_adducts_len *
                    (1 * (topquant_cor))

                  if (good_adducts_len > 0) {
                    chemical_score <- sum(chemical_score * (as.numeric(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts), 2])))
                  }
                }
              }

              names(chemical_score) <- chemicalid[1]

              if (chemical_score > temp_best_score) {
                temp_best_data <- mchemicaldata
                temp_best_score <- chemical_score
              }
            }

            mchemicaldata <- temp_best_data
            chemical_score <- temp_best_score
          } else if (nrow(mchemicaldata) > 1) {
            diff_rt <- abs(min(as.numeric(mchemicaldata$time)) - max(as.numeric(mchemicaldata$time)))

            d1 <- density(mchemicaldata$time, bw = max_diff_rt, from = min(mchemicaldata$time) - 0.001, to = (0.01 + max(mchemicaldata$time)), na.rm = TRUE)
            s1 <- summary(d1$x)

            min_val <- s1[2]
            max_val <- s1[5]

            dup_add <- which(duplicated(mchemicaldata$Adduct) == TRUE)
            if (length(dup_add) > 0) {
              dup_data <- mchemicaldata[dup_add,]
            }

            if (length(which(mchemicaldata$time >= min_val & mchemicaldata$time <= max_val)) > 1) {
              mchemicaldata <- mchemicaldata[which(mchemicaldata$time >= min_val & mchemicaldata$time <= max_val),]

              cur_adducts_with_isotopes <- mchemicaldata$Adduct
              cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")

              good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
              chemical_score <- length(unique(cur_adducts)) *
                good_adducts_len *
                (1 * (topquant_cor))

              if (good_adducts_len > 0) {
                chemical_score <- sum(chemical_score * (10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))))
                chemical_score <- chemical_score[1]
              }

              cur_adducts_with_isotopes <- mchemicaldata$Adduct
              check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])")

              if (length(which(check2 > 0)) > 0) {
                chemical_score <- 100 * chemical_score
              }
            } else {
              cur_adducts_with_isotopes <- mchemicaldata$Adduct
              cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")

              if (diff_rt > max_diff_rt) {
                good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
                chemical_score <- length(unique(cur_adducts)) *
                  good_adducts_len *
                  (1 * (topquant_cor))

                if (good_adducts_len > 0) {
                  chemical_score <- sum(chemical_score * (as.numeric(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts), 2])))
                }
              }
            }
          }

          names(chemical_score) <- chemicalid[1]
          mchemicaldata <- mchemicaldata #[,c(1:11)]
          mchemicaldata <- na.omit(mchemicaldata)
        }
      }

      else {
        cur_adducts_with_isotopes <- mchemicaldata$Adduct
        good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))

        if (good_adducts_len > 0) {
          max_adduct_weight <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))[1]
          chemical_score <- ((10^max_adduct_weight))
          chemical_score <- chemical_score[1]
          good_adduct_index <- which(adduct_weights[, 2] == max_adduct_weight)
          chemical_score <- chemical_score[1]
          mchemicaldata <- mchemicaldata[which(cur_adducts_with_isotopes %in% adduct_weights[good_adduct_index, 1]),]
        } else {
          chemical_score <- 0
          mchemicaldata <- mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust == top_mod[i]),]
          mchemicaldata <- mchemicaldata[order(mchemicaldata$mz),]
        }
      }
    }

    cur_adducts_with_isotopes <- mchemicaldata$Adduct
    cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")

    check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")

    for (a1 in seq_along(check2)) {
      strlength <- attr(check2[[a1]], "match.length")

      if (strlength[1] > -1) {
        count_abundant_form <- length(which(cur_adducts %in% cur_adducts[a1]))

        if (count_abundant_form < 2) {
          mchemicaldata <- mchemicaldata[-a1,]
        }
      }
    }

    if (nrow(mchemicaldata) > 0) {
      conf_level <- get_confidence_stage2(curdata = mchemicaldata, adduct_weights = adduct_weights)
      conf_level <- as.numeric(as.character(conf_level))
    } else {
      conf_level <- 0
    }

    if (nrow(mchemicaldata) > 1) {
      diff_rt <- max(mchemicaldata$time) - min(mchemicaldata$time)

      k_power <- if (diff_rt > max_diff_rt) 10 else 1
      chemical_score <- chemical_score * (1 / ((diff_rt * 0.1) + 1)^k_power)
    } else {
      chemical_score <- 0
    }

    min_chemical_score <- 200 * corthresh * (1 / ((max_diff_rt * 0.1) + 1)^3)

    if (chemical_score > min_chemical_score) {
      chemical_score <- chemical_score * (conf_level^conf_level)
    } else {
      chemical_score <- 0
    }

    if (length(dup_add) > 0) {
      mchemicaldata <- rbind(mchemicaldata, dup_data)
    }

    if (is.na(conf_level) == TRUE) {
      conf_level <- 0
    }

    if (is.na(chemical_score) == TRUE) {
      chemical_score <- 0
    }

    if (chemical_score > best_chemical_score & conf_level > 0) {
      best_chemical_score <- chemical_score
      best_data <- mchemicaldata
    } else if (chemical_score == best_chemical_score) {
      best_chemical_score <- chemical_score
      best_data <- rbind(best_data, mchemicaldata)
    }
  }

  if (best_chemical_score > 0) {
    chemical_score <- best_chemical_score
    mchemicaldata <- best_data
    names(chemical_score) <- chemicalid[1]
  } else {
    chemical_score <- 0
  }
  #######add code for only correlation criteria here


  if (chemical_score <= 1) {
    mchemicaldata <- mchemicaldata_orig
    cur_adducts_with_isotopes <- mchemicaldata$Adduct
    good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))

    if (good_adducts_len > 0) {
      max_adduct_weight <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))
      chemical_score <- 10^max_adduct_weight - 1

      good_adduct_index <- which(adduct_weights[, 2] == max_adduct_weight)
      mchemicaldata <- mchemicaldata[which(cur_adducts_with_isotopes %in% adduct_weights[good_adduct_index, 1]),]
    }
  }

  if (nrow(mchemicaldata) > 0) {
    mchemicaldata <- unique(mchemicaldata)
  } else {
    chemical_score <- 0
  }

  return(list("chemical_score" = chemical_score, "filtdata" = mchemicaldata))
}
