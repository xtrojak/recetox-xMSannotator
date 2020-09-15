compute_confidence_scores <- function(annotation, expected_adducts) {
  #monoisotopic_adduct <- tibble(adduct = "M", charge = 1, mass = 0, molecules = 1)
  # adduct_table <- bind_rows(adduct_table, monoisotopic_adduct)

  # group_by(metabolite) %>%
  # mutate(confidence = if (n() > 1 & all(score < 10) & n_distinct(adducts_with_isotopes) < 2) 0 else confidence) %>%
  # ungroup()

  annotation %>%
    mutate(expected = adduct %in% {{ expected_adducts }}) %>%
    group_by(metabolite) %>%
    mutate(confidence = ifelse(any(expected), 2, 0)) %>%
    ungroup()
}

boost_scores <- function (annotation, boost, mz_tolerance, rt_tolerance) {
  any_near <- function (x, y, tol) any(near(x, y, tol = tol))
  any_rel_near <- function (x, y, tol) any(near(x, y, tol = x * tol))

  boost_by_id <- "metabolite" %in% colnames(boost)
  boost_by_mz <- "mz" %in% colnames(boost)
  boost_by_rt <- "rt" %in% colnames(boost)

  mask <- logical(nrow(annotation))

  if (boost_by_id)
    mask <- annotation$metabolite %in% boost$metabolite
  if (boost_by_rt)
     mask <- mask & sapply(annotation$rt, any_near, boost$rt, rt_tolerance)
  if (boost_by_mz)
    mask <- mask & sapply(annotation$mz, any_rel_near, boost$mz, mz_tolerance)

  annotation$confidence[mask] <- 4
  annotation$score[mask] <- 100 * annotation$score[mask]
  annotation
}

step4 <- function(annotation, expected_adducts, boost_metabolites, mz_tolerance, rt_tolerance) {
  annotation %>%
    compute_confidence_scores(expected_adducts) %>%
    boost_scores(boost_metabolites, mz_tolerance, rt_tolerance) %>%
    mutate(multiple_match = is_nonuinque_mz(mz)) %>%
    print_confidence_distribution()
}


#mutate(
#  expected = adduct %in% expected_adducts,
#  delta_rt = abs(diff(range(rt))),
#  confidence = case_when(
#    any(score < 0.1) ~ 0,
#    !any(adduct %in% adduct_table) ~ 0,
#    # n_distinct(adduct_with_isotopes) < 2 & !any(expected) & any(score < 10) ~ 0,
#    # n_distinct(adduct_with_isotopes) < 2 & any(expected) ~ 1,
#    n_distinct(adduct_with_isotopes) < 2 & !any(expected) ~ 0,
#    n_distinct(adduct_with_isotopes) < 2 & any(expected) ~ 2,
#
#    delta_rt > rt_tolerance & all(is.na(adduct_weight)) ~ 0,
#    delta_rt > rt_tolerance & !all(is.na(adduct_weight)) & !any(expected) &  ~ NA,
#    delta_rt > rt_tolerance & !all(is.na(adduct_weight)) & any(expected) &  ~ NA,
#
#    n_distinct(peak_cluster) > 1 & any(score < 10) ~ 0,
#    n_distinct(peak_cluster) > 1 & !any(score < 10) & any(expected) ~ 2,
#
#    n() < 2 & !any(!is.na(adduct_weight)) ~ 1,
#    n() < 2 & any(!is.na(adduct_weight)) & !any(score < 10) & any(expected) ~ 2,
#    n() < 2 & any(!is.na(adduct_weight)) & !(!any(score < 10) & any(expected)) ~ 0,
#
#    !(n() < 2) & all(molecules > 1) ~ 1,
#    !(n() < 2) & !all(molecules > 1) & ~ # TODO
#  )
#)

#df %>% group_by(metabolite) %>% mutate({
#  curdata <- data.frame()
#  get_confidence_stage4(
#    curdata = data.frame(),
#    max_diff_rt = {{ rt_tolerance }},
#    adduct_weights = data.frame(),
#    filter.by = as.character({{ expected_adducts }}),
#    min_ions_perchem = as.integer({{ min_isp }})
#  )
#})