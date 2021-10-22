#' @export
replace_with_module <- function(module_rt_clust) {
  gsub(as.character(module_rt_clust), pattern = "_[0-9]*", replacement = "")
}

load_adduct_weights <- function() {
  data(adduct_weights)

  adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
  adduct_weights1[1, ] <- c("M+H", 1)
  adduct_weights1[2, ] <- c("M-H", 1)
  adduct_weights <- as.data.frame(adduct_weights1)
  colnames(adduct_weights) <- c("Adduct", "Weight")
  return(adduct_weights)
}

#' @export
extract_possible_isotopes_from_module <- function(...,
                                                  mass_defect_window,
                                                  isop_res_md) {
  query <- tibble(...)

  # Get Module ID
  module_rt_group <- replace_with_module(query$Module_RTclust)
  query_md <- query$mz - round(query$mz)

  md_diff_smaller_window <- abs(isop_res_md$MD - query_md) < mass_defect_window
  is_in_module <- isop_res_md$Module_RTclust == module_rt_group
  extract_indices <- which(md_diff_smaller_window & is_in_module)

  isop_res_md[extract_indices, ]
}

compute_filtered_peak_table <- function(isop_res_md, curmchemdata, mass_defect_window) {
  isop_res_md$Module_RTclust <- replace_with_module(isop_res_md$Module_RTclust)

  isp_masses_mz_data <- purrr::pmap_dfr(
    curmchemdata,
    ~ extract_possible_isotopes_from_module(
      ...,
      mass_defect_window = mass_defect_window,
      isop_res_md = isop_res_md
    )
  )

  colnames(isp_masses_mz_data)[5] <- "AvgIntensity"
  return(isp_masses_mz_data)
}

#' @import plyr
#' @export
compute_chemscore <- function(chemid,
                              mchemdata,
                              isop_res_md,
                              mass_defect_window,
                              corthresh,
                              global_cor,
                              max_diff_rt,
                              adduct_weights,
                              filter.by,
                              max_isp,
                              MplusH.abundance.ratio.check,
                              mass_defect_mode,
                              outloc) {
  chemscoremat <- {}
  curmchemdata <- mchemdata[which(mchemdata$chemical_ID == chemid), ]

  isp_masses_mz_data <- compute_filtered_peak_table(isop_res_md, curmchemdata, mass_defect_window)

  if (is.na(mass_defect_mode) == TRUE) {
    mass_defect_mode <- "pos"
  }

  chem_score <- get_chemscorev1.6.71(
    chemicalid = chemid,
    mchemicaldata = curmchemdata,
    corthresh = corthresh,
    global_cor = global_cor,
    max_diff_rt = max_diff_rt,
    level_module_isop_annot = isp_masses_mz_data,
    adduct_weights = adduct_weights,
    filter.by = filter.by,
    max_isp = max_isp,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
    mass_defect_window = mass_defect_window,
    mass_defect_mode = mass_defect_mode,
    outlocorig = outloc
  )

  if (chem_score$chemical_score >= (-100)) {
    chem_score$filtdata <- chem_score$filtdata[order(chem_score$filtdata$mz), ]
    cur_chem_score <- chem_score$chemical_score
    chemscoremat <- cbind(cur_chem_score, chem_score$filtdata)
    chemscoremat <- as.data.frame(na.omit(chemscoremat))
  }
  return(chemscoremat)
}

#' @import plyr
#' @export
multilevelannotationstep2 <- function(outloc,
                                      max_diff_rt,
                                      mchemdata,
                                      mass_defect_window,
                                      corthresh,
                                      global_cor,
                                      adduct_weights,
                                      max_isp,
                                      MplusH.abundance.ratio.check,
                                      mass_defect_mode,
                                      chemids,
                                      isop_res_md,
                                      filter.by) {

  if (anyNA(adduct_weights)) {
    adduct_weights <- load_adduct_weights()
  }

  chem_score <- lapply(
    chemids,
    compute_chemscore,
    mchemdata = mchemdata,
    isop_res_md = isop_res_md,
    mass_defect_window = mass_defect_window,
    corthresh = corthresh,
    global_cor = global_cor,
    max_diff_rt = max_diff_rt,
    adduct_weights = adduct_weights,
    filter.by = filter.by,
    max_isp = max_isp,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
    mass_defect_mode = mass_defect_mode,
    outloc = outloc
  )

  curchemscoremat <- plyr::ldply(chem_score, rbind)
  return(curchemscoremat)
}
