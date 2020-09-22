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

#compute_scores <- function(df, correlation_threshold, correlation_matrix, adduct_weights, expected_adducts, max_isp, rt_tolerance) {
#  data.frame(
#    mz = mz,
#    time = rt,
#    MatchCategory = ifelse(multiple_match, "Multiple", "Unique"),
#    theoretical.mz = expected_mass,
#    chemical_ID = metabolite,
#    Name = name,
#    Formula = formula,
#    MonoisotopicMass = monoisotopic_mass,
#    Adduct = adduct,
#    ISgroup =,
#    Module_RTclust =,
#    time.y =,
#    mean_int_vec =,
#    MD =
#  )
#
#  as_tibble(df) %>%
#    group_by(metabolite) %>%
#    summarise(
#      data = get_chemscorev1.6.71(
#        chemicalid = metabolite[[1]],
#        mchemicaldata = ,
#        corthresh = {{ correlation_threshold }},
#        global_cor = {{ correlation_matrix }},
#        mzid = paste(mz, rt, sep = "_"),
#        max_diff_rt = rt_tolerance,
#        level_module_isop_annot = ,
#        adduct_weights = as.data.frame({{ adduct_weights }}),
#        filter.by = as.character({{ expected_adducts }}),
#        max_isp = max_isp,
#        MplusH.abundance.ratio.check = ,
#        mass_defect_window = ,
#        mass_defect_mode = "pos",
#      )
#    ) %>%
#    ungroup()
#}