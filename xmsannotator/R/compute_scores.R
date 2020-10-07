compute_isotopes <- function(peaks, query, rt_tolerance, mass_defect_tolerance = 0, abudance_ratio, max_isp) {
  peaks %>%
    dplyr::mutate(
      compound = query$compound,
      isotope_number = round(mz - query$monoisotopic_mass),
      adduct = sprintf("%s_[%+d]", query$adduct, isotope_number),
      formula = sprintf("%s_[%+d]", query$formula, isotope_number),
    ) %>%
    dplyr::filter(
      cluster == query$cluster,
      avg_intensity / query$avg_intensity <= abudance_ratio,
      dplyr::near(rt, query$rt, rt_tolerance),
      dplyr::near(mass_defect, query$mass_defect, mass_defect_tolerance),
      dplyr::between(abs(isotope_number), 1, max_isp)
    ) %>%
    dplyr::select(-isotope_number)
}

compute_scores <- function(annotation, adduct_weights, time_tolerance) {
  annotation <- annotation %>%
    dplyr::filter(is_valid_adduct(adduct, formula)) %>%
    dplyr::left_join(adduct_weights, by = "adduct")
  isotopes <- annotation %>%
    dplyr::filter(!is.na(adduct_weight)) %>%
    purrr::pmap_dfr(~tibble()) %>%
    dplyr::distinct()
  annotation <- dplyr::bind_rows(annotation, isotopes) %>%
    dplyr::group_by(compound) %>%
    dplyr::mutate(score = 10^max(-Inf, adduct_weight, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(multiple_match = is_nonuinque_mass(mz))
}
