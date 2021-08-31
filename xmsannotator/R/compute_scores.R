#' Assign abundance ratio of the second most abundant isotope/isotopologue for each entry in a table
#'
#' @param isotopes A table with annotated peaks filtered by adducts
#' @return The input table with extra column containing the abundance ratio of the second most abundant isotope/isotopologue
#'
#' @import dplyr
#' @import purrr
assign_isotope_abundances <- function(isotopes) {
  molecules <- distinct(isotopes, molecular_formula)
  molecules <- mutate(molecules, abundance_ratio = map_dbl(molecular_formula, compute_abundance_ratio))
  isotopes <- inner_join(isotopes, molecules, by = "molecular_formula")
}

#' Get abundance ratio of the second most abundant isotope/isotopologue
#'
#' @param formula A string containing molecular or empirical formula of a compound
#'
#' @return Normalized abundance ratio of the second most abundant isotope/isotopologue
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

#' Match isotopes from peak table to annotated peaks
#'
#' @param ... A single isotope passed as a list of columns
#' @param intensity_deviation_tolerance A number. A threshold by which an intensity ratio of isotope
#' @param peaks A peak table - `isp_masses_mz_data` from the original tool's multilevelannotation step 1
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound
#' @param max_isp A number. Maximal number of unique isotopes of a single compound
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance
#'  to a molecular peak may exceed its relative isotopic abundance. The default value was hardcoded by the
#'  original author.
#'
#' @return A table of matching isotopes from `peaks` table
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

#' Compute a chemical (confidence) score for each annotation
#'
#' @param annotation A table with annotated peaks
#' @param adduct_weights A weight-by-adduct table
#' @param intensity_deviation_tolerance A number. A threshold by which an intensity ratio of isotope
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound
#' @param max_isp A number. Maximal number of unique isotopes of a single compound
#' @param peaks A feature table containing mz, rt, rt cluster, and mean intensity of each peak
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
