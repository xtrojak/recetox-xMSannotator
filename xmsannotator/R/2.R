make_isotopes <- function(peaks, query, rt_tolerance, mass_defect_tolerance, abudance_ratio, max_isp) {
  peaks %>%
    mutate(
      metabolite = query$metabolite,
      isotope_number = round(mz - query$monoisotopic_mass),
      adduct = sprintf("%s_[%+d]", query$adduct, isotope_number),
      formula = sprintf("%s_[%+d]", query$formula, isotope_number),
    ) %>%
    filter(
      cluster == query$cluster,
      avg_intensity / query$avg_intensity <= abudance_ratio,
      near(rt, query$rt, rt_tolerance),
      near(mass_defect, query$mass_defect, mass_defect_tolerance),
      between(abs(isotope_number), 1, max_isp)
    ) %>%
    select(-isotope_number)
}

compute_scores <- function(annotation, adduct_weights, rt_tolerance, mass_defect_tolerance) {
  annotation <- annotation %>%
    filter(is_valid_adduct(adduct, formula)) %>%
    left_join(adduct_weights, by = "adduct")
  isotopes <- annotation %>%
    filter(!is.na(weight)) %>%
    pmap_dfr(~tibble()) %>%
    distinct()

  bind_rows(annotation, isotopes) %>%
    group_by(metabolite) %>%
    mutate(score = 10^max(-Inf, weight, na.rm = TRUE)) %>%
    ungroup()
}