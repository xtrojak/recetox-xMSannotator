mass_defect <- function (x, precision) {
  cut(x %% 1, seq(0, 1, precision), labels = FALSE)
}

print_confidence_distribution <- function(annotation) {
  confidence_distribution_across_compounds <- dplyr::filter(annotation, !duplicated(compound))
  confidence_distribution_across_compounds <- dplyr::count(confidence_distribution_across_compounds, confidence)

  confidence_distribution_across_formulas <- dplyr::filter(annotation, !duplicated(molecular_formula))
  confidence_distribution_across_formulas <- dplyr::count(confidence_distribution_across_formulas, confidence)

  print("Confidence level distribution for unique compounds")
  print(confidence_distribution_across_compounds)
  print("Confidence level distribution for unique formulas")
  print(confidence_distribution_across_formulas)

  invisible(annotation)
}

advanced_annotation <- function(peak_table, compound_table, adduct_table,
  adduct_weights = tibble::tibble(adduct = c("M+H", "M-H"), adduct_weight = c(5, 5)),
  mass_tolerance = 5e-6,
  time_tolerance = 10,
  correlation_threshold = 0.7,
  deep_split = 2,
  min_cluster_size = 10,
  network_type = "unsigned",
  expected_adducts = character(),
  boost_compounds = tibble::tibble(compound = character(), mz = numeric(), rt = numeric()),
  redundancy_filtering = TRUE,
  n_workers = parallel::detectCores()
) {
  if (is.numeric(n_workers) && n_workers > 1) {
    WGCNA::allowWGCNAThreads(n_workers)
  }

  peak_table <- dplyr::select(peak_table, peak, mz, rt, dplyr::starts_with("intensity"))
  peak_table <- dplyr::distinct(peak_table, peak, .keep_all = TRUE)
  peak_table <- tidyr::pack(peak_table, intensity = dplyr::starts_with("intensity"))

  intensity_matrix <- t(dplyr::pull(peak_table, intensity))
  intensity_matrix <- magrittr::set_colnames(intensity_matrix, dplyr::pull(peak_table, peak))
  correlation_matrix <- WGCNA::cor(intensity_matrix, use = "p", method = "p")

  clustering <- peak_clustering(
    intensity_matrix = intensity_matrix,
    correlation_matrix = correlation_matrix,
    correlation_threshold = correlation_threshold,
    min_cluster_size = min_cluster_size,
    network_type = network_type,
    deep_split = deep_split
  )

  peak_info <- dplyr::transmute(
    peak_table,
    peak,
    cluster = clustering,
    intensity = rowMeans(intensity, na.rm = TRUE),
    mass_defect = mass_defect(mz, precision = 0.01)
  )

  annotation <- simple_annotation(peak_table, adduct_table, compound_table, mass_tolerance)
  annotation <- dplyr::inner_join(annotation, peak_info, by = "peak")

  annotation <- compute_scores(annotation, adduct_weights, time_tolerance)
  # annotation <- evaluate_pathways(annotation, pathways, score_threshold = 0.1)
  annotation <- compute_confidence_levels(annotation, expected_adducts)
  annotation <- boost_scores(annotation, boost_compounds, mass_tolerance, time_tolerance)
  print_confidence_distribution(annotation)

  if (redundancy_filtering) {
    annotation <- remove_duplicates(annotation, score_threshold = 0)
    print_confidence_distribution(annotation)
  }

  return(annotation)
}
