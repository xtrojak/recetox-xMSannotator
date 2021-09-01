#' Assign the abundance ratio of the second most abundant isotope/isotopologue for each molecule in a table
#'
#' @param isotopes A table with annotated peaks
#'
#' @return The input table with an extra column containing the normalized abundance ratio of the second most abundant isotope/isotopologue
#'
#' @import dplyr
#' @import purrr
assign_isotope_abundances <- function(isotopes) {
  molecules <- distinct(isotopes, molecular_formula)
  molecules <- mutate(molecules, abundance_ratio = map_dbl(molecular_formula, compute_abundance_ratio))
  isotopes <- inner_join(isotopes, molecules, by = "molecular_formula")
}

#' Compute abundance ratio of the second most abundant isotope/isotopologue
#'
#' @param formula A string containing molecular or empirical formula of a compound
#'
#' @return Normalized abundance ratio of the second most abundant isotope/isotopologue. The normalization is performed
#'  such that the most abundant isotope has a unit abundance and abundances of the rest of the isotopes are
#'  represented as a share of this unit abundance
#'
#' @import dplyr
#' @importFrom rcdk get.formula get.isotopes.pattern
compute_abundance_ratio <- function(formula) {
  formula <- get.formula(formula)
  isotopes <- get.isotopes.pattern(formula, minAbund = 0.001)
  isotopes <- as_tibble(isotopes)
  isotopes <- arrange(isotopes, desc(abund))
  sec_most_abundant <- isotopes$abund[2]
  sec_most_abundant <- max(c(sec_most_abundant, 0), na.rm = TRUE)
}

#' Match isotopes from peak table to annotated peaks based on rt cluster, rt difference, peak intensity, and mass defect.
#'
#' @param ... A single isotope passed as a list of columns. See [purrr::map()] for more details
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio
#' @param peaks A peak table containing a peak identifier (unique number), mean intensity, module, and rt cluster
#'  of each identified peak
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound
#' @param max_isp Maximal number of unique isotopes of a single compound
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance
#'
#' @return A table of identified isotopes with an extra column: `isotopic deviation`.
#'  The extra column represents a nominal mass (mass number) difference between a given isotope and the most abundant
#'  isotope of the same molecule
#'
#' @import dplyr
#' @importFrom rlang .data
compute_isotopes <- function(...,
                             intensity_deviation_tolerance = 0.1,
                             peaks,
                             mass_defect_tolerance = 0,
                             max_isp,
                             rt_tolerance) {
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
    RTclust == query$RTclust,
    mean_intensity / query$mean_intensity <= query$abundance_ratio + intensity_deviation_tolerance,
    near(rt, query$rt, rt_tolerance),
    near(mass_defect, query$mass_defect, mass_defect_tolerance),
    between(abs(isotopic_deviation), 1, max_isp)
  )
}

#' Compute a confidence score for each annotation
#'
#' @param annotation A table with annotated peaks
#' @param adduct_weights A weight-by-adduct table
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound
#' @param max_isp Maximal number of unique isotopes of a single compound
#' @param peaks A peak table containing a peak identifier (unique number), mean intensity, module, and rt cluster
#'  of each identified peak
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance
#'
#' @import dplyr
#' @importFrom rlang .data
compute_scores <- function(annotation,
                           adduct_weights,
                           intensity_deviation_tolerance = 0.1,
                           mass_defect_tolerance,
                           max_isp,
                           peaks,
                           rt_tolerance) {
  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))

  isotopes <- semi_join(annotation, adduct_weights, by = "adduct")
  isotopes <- assign_isotope_abundances(isotopes)
  # This can be parallelized on `group_split(group_by(isotopes, molecular_formula))`
  isotopes <- purrr::pmap_dfr(isotopes,
                              ~compute_isotopes(..., peaks = peaks,
                                                rt_tolerance = rt_tolerance,
                                                intensity_deviation_tolerance = intensity_deviation_tolerance,
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
