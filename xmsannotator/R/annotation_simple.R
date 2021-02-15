is_nonuinque_mass <- function(mass) {
  freq <- tibble::as_tibble_col(mass)
  freq <- dplyr::add_count(freq, value)
  freq <- dplyr::transmute(freq, n > 1)
  dplyr::pull(freq)
}

forms_valid_adduct_pair <- function(molecular_formula, adduct_formula) {
  has_oxygen <- check_element(molecular_formula, "O") > 0
  is_water_adducts <- adduct_formula %in% c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
  !is_water_adducts | has_oxygen
}

simple_annotation <- function(peak_table, adduct_table, compound_table, mass_tolerance = 5e-6) {
  peak_table <- dplyr::select(peak_table, peak, mz, rt)
  peak_table <- dplyr::distinct(peak_table, peak, .keep_all = TRUE)

  adduct_table <- dplyr::select(adduct_table, adduct, charge, mass, n_molecules)
  adduct_table <- dplyr::distinct(adduct_table, adduct, .keep_all = TRUE)

  compound_table <- dplyr::select(compound_table, compound, monoisotopic_mass, molecular_formula)
  compound_table <- dplyr::distinct(compound_table, compound, .keep_all = TRUE)

  annotation <- annotate_by_mass(peak_table, adduct_table, compound_table, mass_tolerance)
  annotation <- dplyr::left_join(annotation, peak_table, by = "peak")
  annotation <- dplyr::left_join(annotation, compound_table, by = "compound")

  # NOTE: the following line is commented out because it only checks compound's formula instead of compound-adduct combination
  annotation <- dplyr::filter(annotation, check_golden_rules(molecular_formula))
  annotation <- dplyr::filter(annotation, forms_valid_adduct_pair(molecular_formula, adduct))
  annotation <- dplyr::mutate(annotation, multiple_match = is_nonuinque_mass(mz))

  return(annotation)
}
