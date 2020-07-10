multilevelannotationstep2 <- function(outloc1, list_number, adduct_table, adduct_weights) {
    setwd(outloc1)

    load("step1_results.Rda")
    load("global_cor.Rda")
    unlink("allmatches_with_isotopes.txt")

    outloc <- outloc1
    if (is.na(max.rt.diff) == FALSE) {
      max_diff_rt <- max.rt.diff
    }

    if (list_number > length(chemids_split)) {
      return(0)
    }

    if (list_number > num_sets) {
      return(0)
    }

    outloc1 <- paste0(outloc, "/stage2/")
    suppressWarnings(dir.create(outloc1))
    setwd(outloc1)

    cnames <- colnames(mchemdata)
    cnames[2] <- "time"
    colnames(mchemdata) <- as.character(cnames)
    mchemdata$mz <- as.numeric(as.character(mchemdata$mz))
    mchemdata$time <- as.numeric(as.character(mchemdata$time))


    chem_score <- lapply(chemids_split[[list_number]], function(j) {
      chemid <- chemids[j]
      chemscoremat <- { }
      curmchemdata <- mchemdata[which(mchemdata$chemical_ID == chemid),]
      curmchemdata$mz <- as.numeric(as.character(curmchemdata$mz))
      curmchemdata$time <- as.numeric(as.character(curmchemdata$time))
      curmchemdata <- as.data.frame(curmchemdata)
      curmchemdata$Module_RTclust <- gsub(curmchemdata$Module_RTclust, pattern = "_[0-9]*", replacement = "")
      isop_res_md$Module_RTclust <- gsub(isop_res_md$Module_RTclust, pattern = "_[0-9]*", replacement = "")
      isp_masses_mz_data <- lapply(seq_along(curmchemdata$mz), function(m)
        {
        isp_group <- as.character(curmchemdata$ISgroup[m])
        module_rt_group <- as.character(curmchemdata$Module_RTclust[m])
        module_rt_group <- gsub(module_rt_group, pattern = "_[0-9]*", replacement = "")
        query_md <- curmchemdata$mz[m] - round(curmchemdata$mz[m])
        put_isp_masses_curmz_data <- isop_res_md[which(abs((isop_res_md$MD) - (query_md)) < mass_defect_window & isop_res_md$Module_RTclust == module_rt_group),]
        put_isp_masses_curmz_data <- as.data.frame(put_isp_masses_curmz_data)
        return(put_isp_masses_curmz_data)
      })
      isp_masses_mz_data <- ldply(isp_masses_mz_data, rbind)
      cnames <- colnames(isp_masses_mz_data)
      cnames[5] <- "AvgIntensity"
      colnames(isp_masses_mz_data) <- cnames
      isp_masses_mz_data <- as.data.frame(isp_masses_mz_data)
      if (is.na(mass_defect_mode) == TRUE) {
        mass_defect_mode <- "pos"
      }

      isp_masses_mz_data$mz <- as.numeric(as.character(isp_masses_mz_data$mz))
      isp_masses_mz_data$time <- as.numeric(as.character(isp_masses_mz_data$time))
      isp_masses_mz_data <- as.data.frame(isp_masses_mz_data)

      chem_score <- get_chemscorev1.6.71(chemicalid = chemid, mchemicaldata = curmchemdata, corthresh = corthresh, global_cor = global_cor, mzid, max_diff_rt = max_diff_rt,
                                         level_module_isop_annot = isp_masses_mz_data, adduct_table = adduct_table, adduct_weights = adduct_weights, filter.by, max_isp = max_isp,
                                         MplusH.abundance.ratio.check = MplusH.abundance.ratio.check, mass_defect_window = mass_defect_window, mass_defect_mode = mass_defect_mode, outlocorig = outloc)


      if (length(chem_score) > 0) {
        if (chem_score$chemical_score >= (-100)) {
          chem_score$filtdata <- chem_score$filtdata[order(chem_score$filtdata$mz),]

          if (nrow(chem_score$filtdata) > 0) {
            cur_chem_score <- chem_score$chemical_score
            chemscoremat <- cbind(cur_chem_score, chem_score$filtdata)
          }

          chemscoremat <- na.omit(chemscoremat)
          chemscoremat <- as.data.frame(chemscoremat)
        }
      }

      return(chemscoremat)
    })

    chem_score2 <- chem_score[which(chem_score != "NULL")]
    curchemscoremat <- plyr::ldply(chem_score2, rbind)

    cur_fname <- paste0("chem_score", list_number, ".Rda")
    save(curchemscoremat, file = cur_fname)
  }
