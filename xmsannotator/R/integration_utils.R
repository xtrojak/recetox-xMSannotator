construct_adduct_formula <- function(formula, mass_number_difference) {
  suffix <- sprintf("[%+1.0f]", mass_number_difference)
  formula <- paste(formula, suffix, sep = "_")
  return(formula)
}

#' Reformat annotation table with isotopes from `main` to be compatible with `master's` chemical score computing.
#'
#' @import dplyr
#' @importFrom magittr %>%
#' @importFrom rlang .data
#' @export
reformat_annotation_table <- function(annotation, supplementary_data) {
  master_annotation <- tibble(.rows = nrow(annotation))

  master_annotation <- master_annotation %>%
    mutate(
      Module_RTclust = paste(annotation$module, annotation$rt_cluster, sep = "_"),
      mz = annotation$mz,
      time = annotation$rt,
      Formula = annotation$molecular_formula,
      Adduct = if_else(annotation$mass_number_difference == 0, annotation$adduct, construct_adduct_formula(annotation$adduct, annotation$mass_number_difference)),
      time.y = annotation$rt,
      ISgroup = character(nrow(annotation)), # dummy column due to high complexity of reproducing the values
      mean_int_vec = annotation$mean_intensity,
      MD = as.numeric(paste(0, annotation$mass_defect, sep = "."))
    )

  return(master_annotation)
}