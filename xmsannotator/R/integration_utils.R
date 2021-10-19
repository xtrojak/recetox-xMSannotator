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
reformat_annotation_table <- function(annotation) {
  master_annotation <- tibble(.rows = nrow(annotation))

  master_annotation <- master_annotation %>%
    mutate(
      Module_RTclust = paste(annotation$module, annotation$rt_cluster, sep = "_"),
      mz = annotation$mz,
      time = annotation$rt,
      MatchCategory = case_when(
        is.na(annotation$multiple_match) ~ "-",
        isTRUE(annotation$multiple_match) ~ "Multiple",
        !isTRUE(annotation$multiple_match) ~ "Unique"
      ),
      theoretical.mz = annotation$expected_mass,
      chemical_ID = paste("Formula", annotation$compound, sep = "_"),
      Name = "-",
      Formula = annotation$molecular_formula,
      MonoisotopicMass = annotation$monoisotopic_mass,
      Adduct = if_else(annotation$mass_number_difference == 0, annotation$adduct, construct_adduct_formula(annotation$adduct, annotation$mass_number_difference)),
      ISgroup = "-", # dummy column due to high complexity of reproducing the values
      time.y = annotation$rt,
      mean_int_vec = annotation$mean_intensity,
      MD = as.numeric(sprintf("0.%1.0f", annotation$mass_defect))
    )

  return(master_annotation)
}
