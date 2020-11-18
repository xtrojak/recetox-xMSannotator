is_nonuinque_mass <- function(mass) {
  tibble::as_tibble_col(mass) %>%
    dplyr::add_count(value) %>%
    dplyr::transmute(n > 1) %>%
    dplyr::pull()
}

simple_annotation <- function(peaks, compounds, adducts, mass_tolerance = 5e-6) {
  peaks <- validate_table(
    table = peaks,
    primary_key = c("mz", "rt")
  ) %>%
    dplyr::mutate(peak = lexicographic_rank(mz, rt)) %>%
    dplyr::distinct(peak, .keep_all = TRUE)

  adducts <- validate_table(
    table = adducts,
    primary_key = "adduct",
    expected_columns = c("charge", "mass", "molecules")
  ) %>%
    dplyr::rename(mode = charge, multiplicity = molecules)

  compounds <- validate_table(
    table = compounds,
    primary_key = "compound",
    expected_columns = c("monoisotopic_mass", "formula")
  )

  annotation <- annotate_by_mass(peaks, adducts, compounds, mass_tolerance) %>%
    dplyr::left_join(peaks, by = "peak") %>%
    dplyr::left_join(compounds, by = "compound") %>%
    # dplyr::filter(check_golden_rules(formula)) %>%
    # dplyr::filter(is_valid_adduct(adduct, formula)) %>%
    dplyr::mutate(multiple_match = is_nonuinque_mass(mz))
}
