is_nonuinque_mass <- function(mass) {
  tibble::as_tibble_col(mass) %>%
    dplyr::add_count(value) %>%
    dplyr::transmute(n > 1) %>%
    dplyr::pull()
}

merge_by_mass <- function(peaks, compounds, mass_tolerance) {
  i_indices <- seq_len(nrow(peaks))
  j_indices <- seq_len(nrow(compounds))

  indices <- tidyr::expand_grid(i = i_indices, j = j_indices) %>%
    dplyr::filter(dplyr::near(peaks$mz[i], compounds$expected_mass[j], mass_tolerance))
  gc() # call garbage collector after heavy memory load

  peaks <- dplyr::slice(peaks, indices$i)
  compounds <- dplyr::slice(compounds, indices$j)

  dplyr::bind_cols(peaks, compounds)
}

simple_annotation <- function(peaks, compounds, adducts, mass_tolerance = 5e-6) {
  peaks <- validate_table(
    table = peaks,
    primary_key = c("mz", "rt")
  )

  adducts <- validate_table(
    table = adducts,
    primary_key = "adduct",
    expected_columns = c("charge", "mass", "molecules")
  )

  compounds <- validate_table(
    table = compounds,
    primary_key = "compound",
    expected_columns = c("monoisotopic_mass", "formula")
  )

  compounds <- tidyr::expand_grid(compounds, adducts) %>%
    dplyr::filter(check_golden_rules(formula)) %>%
    dplyr::filter(is_valid_adduct(adduct, formula)) %>%
    dplyr::mutate(expected_mass = (molecules * monoisotopic_mass + mass) / charge) %>%
    dplyr::select(adduct, compound, formula, monoisotopic_mass, expected_mass)
  annotation <- merge_by_mass(peaks, compounds, mass_tolerance) %>%
    dplyr::mutate(multiple_match = is_nonuinque_mass(mz))
}