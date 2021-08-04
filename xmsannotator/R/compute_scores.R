# Get abundance ratio of the second most abundant isotopologue
#' @import Rdisop
compute_abundance_ratio <- function(formula) {
  molecule <- getMolecule(formula)
  isotope <- getIsotope(molecule, 2)
  abundance_ratio <- isotope[2]
}

#' @import dplyr
#' @importFrom rlang .data
compute_isotopes <- function(..., peaks, rt_tolerance, mass_defect_tolerance = 0, max_isp) {
  query <- tibble(...)
  isotopes <- mutate(
    peaks,
    isotopic_deviation = round(mz - query$expected_mass),
    compound = query$compound,
    adduct = query$adduct,
    molecular_formula = query$molecular_formula
  )
  isotopes <- filter(
    isotopes,
    cluster == query$cluster,
    mean_intensity / query$mean_intensity <= query$abundance_ratio,
    near(rt, query$rt, rt_tolerance),
    near(mass_defect, query$mass_defect, mass_defect_tolerance),
    between(abs(isotopic_deviation), 1, max_isp)
  )
}

#' @import dplyr
#' @importFrom rlang .data
compute_scores <- function(annotation, adduct_weights, mass_defect_tolerance, max_isp, peaks, rt_tolerance) {
  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))

  isotopes <- semi_join(annotation, adduct_weights, by = "adduct")
  isotopes <- mutate(isotopes, abundance_ratio = map_dbl(molecular_formula, compute_abundance_ratio))
  isotopes <- purrr::pmap_dfr(isotopes,
                              ~compute_isotopes(..., peaks = peaks,
                                                rt_tolerance = rt_tolerance,
                                                mass_defect_tolerance = mass_defect_tolerance,
                                                max_isp = max_isp))
  isotopes <- distinct(isotopes)

  annotation <- bind_rows(annotation, isotopes)
  annotation <- left_join(annotation, adduct_weights, by = "adduct")
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, score = 10^max(-Inf, weight, na.rm = TRUE))
  annotation <- ungroup(annotation)
  annotation <- compute_multiple_matches(annotation)
  annotation
}
