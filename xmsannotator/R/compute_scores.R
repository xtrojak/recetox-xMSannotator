compute_isotopes <- function(peaks, query, rt_tolerance, mass_defect_tolerance = 0, abudance_ratio, max_isp) {
  isotopes <- dplyr::mutate(
    peaks,
    compound = query$compound,
    isotope_number = round(mz - query$monoisotopic_mass),
    adduct = sprintf("%s_[%+d]", query$adduct, isotope_number),
    formula = sprintf("%s_[%+d]", query$formula, isotope_number),
  )
  isotopes <- dplyr::filter(
    isotopes,
    cluster == query$cluster,
    avg_intensity / query$avg_intensity <= abudance_ratio,
    dplyr::near(rt, query$rt, rt_tolerance),
    dplyr::near(mass_defect, query$mass_defect, mass_defect_tolerance),
    dplyr::between(abs(isotope_number), 1, max_isp)
  )
  isotopes <- dplyr::select(isotopes, -isotope_number)
}

compute_scores <- function(annotation, adduct_weights, time_tolerance) {
  annotation <- dplyr::filter(annotation, forms_valid_adduct_pair(molecular_formula, adduct))
  annotation <- dplyr::left_join(annotation, adduct_weights, by = "adduct")

  isotopes <- dplyr::filter(annotation, !is.na(adduct_weight))
  isotopes <- purrr::pmap_dfr(isotopes, ~tibble::tibble())
  isotopes <- dplyr::distinct(isotopes)

  annotation <- dplyr::bind_rows(annotation, isotopes)
  annotation <- dplyr::group_by(annotation, compound)
  annotation <- dplyr::mutate(annotation, score = 10^max(-Inf, adduct_weight, na.rm = TRUE))
  annotation <- dplyr::ungroup(annotation)
  annotation <- dplyr::mutate(annotation, multiple_match = is_nonuinque_mass(mz))
}
