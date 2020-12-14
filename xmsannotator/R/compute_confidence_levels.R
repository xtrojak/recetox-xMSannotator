compute_confidence_levels <- function(annotation, expected_adducts) {
  # monoisotopic_adduct <- tibble::tibble(adduct = "M", charge = 1, mass = 0, molecules = 1)
  # adduct_table <- dplyr::bind_rows(adduct_table, monoisotopic_adduct)

  # dplyr::group_by(compound) %>%
  # dplyr::mutate(confidence = if (n() > 1 & all(score < 10) & n_distinct(adducts_with_isotopes) < 2) 0 else confidence) %>%
  # dplyr::ungroup()

  annotation %>%
    dplyr::mutate(expected = adduct %in% {{ expected_adducts }}) %>%
    dplyr::group_by(compound) %>%
    dplyr::mutate(confidence = ifelse(any(expected), 2, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(multiple_match = is_nonuinque_mass(mz))
}

boost_scores <- function(annotation, boost_compounds, mass_tolerance, time_tolerance) {
  boost_compounds <- dplyr::select(boost_compounds, dplyr::any_of(c("compound", "mz", "rt")))

  any_near <- function(x, y, tol) any(near(x, y, tol = tol))
  any_rel_near <- function(x, y, tol) any(near(x, y, tol = x * tol))

  boost_by_id <- "compound" %in% colnames(boost_compounds)
  boost_by_mz <- "mz" %in% colnames(boost_compounds)
  boost_by_rt <- "rt" %in% colnames(boost_compounds)

  mask <- logical(nrow(annotation))

  if (boost_by_id)
    mask <- annotation$compound %in% boost_compounds$compound
  if (boost_by_rt)
    mask <- mask & sapply(annotation$rt, any_near, boost_compounds$rt, time_tolerance)
  if (boost_by_mz)
    mask <- mask & sapply(annotation$mz, any_rel_near, boost_compounds$mz, mass_tolerance)

  annotation$confidence[mask] <- 4
  annotation$score[mask] <- 100 * annotation$score[mask]
  annotation
}
