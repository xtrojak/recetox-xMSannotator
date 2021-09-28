#' Compute a confidence score for each annotation.
#'
#' @param annotation A table with annotated peaks.
#' @param adduct_weights A weight-by-adduct table.
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio.
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound.
#' @param peaks A peak table containing a peak identifier (unique number), mean intensity, module, and rt cluster
#'  of each identified peak.
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance.
#'
#' @import dplyr
#' @import purrr
#' @importFrom rlang .data
compute_scores <- function(annotation,
                           adduct_weights,
                           intensity_deviation_tolerance = 0.1,
                           mass_defect_tolerance = 0,
                           peaks,
                           rt_tolerance) {
  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))

  isotopes <- semi_join(annotation, adduct_weights, by = "adduct")
  # This can be parallelized on `group_split(group_by(isotopes, molecular_formula))`
  isotopes <- purrr::pmap_dfr(isotopes,
                              ~compute_isotopes(..., peaks = peaks,
                                                rt_tolerance = rt_tolerance,
                                                intensity_deviation_tolerance = intensity_deviation_tolerance,
                                                mass_defect_tolerance = mass_defect_tolerance))

  annotation <- bind_rows(annotation, isotopes)
  annotation <- left_join(annotation, adduct_weights, by = "adduct")
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, score = 10^max(-Inf, weight, na.rm = TRUE))
  annotation <- ungroup(annotation)
  annotation <- compute_multiple_matches(annotation)
  annotation
}
