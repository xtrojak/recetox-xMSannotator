advanced_annotation <- function(peak_intensity_table, metabolites, workers) {
  WGCNA::allowWGCNAThreads(workers)

  peak_intensity_table <- as_tibble(peak_intensity_table) %>%
    distinct(peak_intensity_table, mz, rt, .keep_all = TRUE)
  metabolites <- dplyr::distinct(metabolites)
}

compute_score <- function(df, adduct_weights) {
  # TODO line 11-16

  isotopes <- Rdisop::getMolecule(formula)$isotopes[[1]]
  abundance_ratio <- isotopes[2, isotopes[2,] >= 0.001]

  # TODO: line 26-41 remove water adduct
  # TODO: line 47-160
  # TODO: ..186 group by rt with isotopes

  dt %>%
    add_count(rt_group, sort = TRUE)
}

compute_scores <- function(df, adduct_weights) {
  df %>%
    group_by(id) %>%
    mutate(score = { }) %>%
    ungroup(id)
}

.pathway_evaluation <- function(dt, metabolites, adduct_weights, score_threshold) {
  key <- function(x) gsub(pattern = "_.*", replacement = "", x)

  dt <- dt %>%
    mutate(adduct_key = key(adduct), formula_key = key(formula)) %>%
    inner_join(metabolites, by = c("adduct_key", "formula_key"))

  hmdb_bad <- c("HMDB29244", "HMDB29245", "HMDB29246")
  dt <- filter(dt, !(chemical_id %in% hmdb_bad))

  # TODO: extract `pathways` from metabolite DB
  semi_join(metabolites, dt, by = pathway)

  filter(dt, score >= score_threshold)
}

print_confidence_distribution <- function(df) {
  unique_chemical_confidence <- df %>%
    select(chemical_id, confidence) %>%
    distinct(chemical_id, .keep_all = TRUE) %>%
    pull(confidence)
  unique_formula_confidence <- df %>%
    select(formula, confidence) %>%
    distinct(formula, .keep_all = TRUE) %>%
    pull(confidence)

  print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
  print(table(unique_chemical_confidence))

  print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
  print(table(unique_formula_confidence))

  invisible(df)
}

step5 <- function(result, adduct_weights, score_threshold = 0) {
  mz_frequency <- table(result$mz)
  nonunique_mzs <- names(mz_frequency)[mz_frequency > 1]

  bad_indices <- purrr::map_int(nonunique_mzs, function(value) {
    subset <- result[result$mz %in% value,]

    good_adduct <- subset$adduct %in% adduct_weights$adduct
    subset$score[good_adduct] <- subset$score[good_adduct] * 100

    indices <- which(result$mz %in% value)
    indices <- indices[subset$score != max(subset$score, na.rm = TRUE)]

    for (chemical_id in subset$chemical_id[indices]) {
      ind <- which(result$chemical_id %in% chemical_id)
      result$score[ind] <- (length(ind) - 1) * result$score[ind[1]] / length(ind)
    }

    indices
  })

  result <- result %>%
    slice(-bad_indices) %>%
    filter(score >= score_threshold) %>%
    mutate(multiple_match = is_nonuinque_mz(mz))


  print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
  # TODO: print(table(unique_chemical_confidence))

  print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
  # TODO: print(table(unique_formula_confidence))

  result
}

compute_multiple_matches <- function (df) {
  mutate(df, multiple_match = is_nonuinque_mz(mz))
}

compute_delta_mz <- function(df) {
  # NOTE: I changed the computation of delta_ppm becaouse I think the previos
  # version behaved incorectly (see lines 130-140 in multilevelannotationstep4).
  mutate(df, delta_ppm = round(10^6 * abs(theoretical_mz - mz) / theoretical_mz, 2))
}

boost_scores <- function(df, boost, mz_diff, rt_diff) {
  df <- rename_with(df, tolower)

  match_peak_by_id <- function(id) id %in% boost$id
  match_peak_by_mz <- function(mz) abs(mz - boost$mz) < mz * mz_diff
  match_peak_by_rt <- function(rt) abs(rt - boost$rt) < rt_diff

  boosted <- if (all(c("id", "mz", "rt") %in% colnames(boost))) {
    match_peak_by_id(id, boost$id)
  } else if ("id" %in% colnames(boost)) {
    mapply(function(id, mz, rt) {
      any(match_peak_by_id(id) & match_peak_by_mz(mz) & match_peak_by_rt(rt))
    }, id = df$id, mz = df$mz, rt = df$rt)
  } else {
    stop("Boost table must atleast contain one column called ID")
  }

  mutate(df,
    score = ifelse(boosted, 100 * score, score),
    confidence = ifelse(boosted, 4, confidence)
  )
}

step4 <- function(df, adduct_table, adduct_weights, filter_by, boost, mz_diff, rt_diff) {
  df %>%
    compute_confidence_scores(
      adduct_table = adduct_table,
      adduct_weights = adduct_weights,
      filter_by = filter_by
    ) %>%
    boost_scores(
      boost = boost,
      mz_diff = mz_diff,
      rt_diff = rt_diff
    ) %>%
    compute_multiple_matches() %>%
    print_confidence_distribution()
}

compute_confidence_scores <- function (df, adduct_table, adduct_weights, filter_by) {
  monoisotopic <- list(
    adduct = "M",
    molecules = 1,
    adduct_charge = 1,
    adduct_mass = 0,
    adduct_mode = NA,
    adduct_type ="S"
  )
  adduct_table <- bind_rows(adduct_table, monoisotopic)

  df %>%
    group_by(chemical_id) %>%
    mutate(confidence = 0) %>%
    ungroup()
}
