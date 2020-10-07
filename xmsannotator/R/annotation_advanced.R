lexicographic_rank <- function (...) {
  .o <- order(...)
  .x <- cbind(...)[.o, ]
  cumsum(!duplicated(.x))[order(.o)]
}

mass_defect <- function (x, precision) {
  cut(x %% 1, seq(0, 1, precision), labels = FALSE)
}

print_confidence_distribution <- function(annotation) {
  confidence_distribution_across_compounds <- annotation %>%
    filter(!duplicated(compound)) %>%
    count(confidence)
  confidence_distribution_across_formulas <- annotation %>%
    filter(!duplicated(formula)) %>%
    count(confidence)

  print("Confidence level distribution for unique compounds")
  print(confidence_distribution_across_compounds)
  print("Confidence level distribution for unique formulas")
  print(confidence_distribution_across_formulas)

  invisible(annotation)
}

advanced_annotation <- function(peaks, compounds, adducts,
  adduct_weights = tibble::tibble(adduct = c("M+H", "M-H"), adduct_weight = c(5, 5)),
  mass_tolerance = 5e-6,
  time_tolerance = 10,
  correlation_threshold = 0.7,
  deep_split = 2,
  min_cluster_size = 10,
  network_type = "unsigned",
  expected_adducts = character(),
  boost = tibble::tibble(compound = character(), mz = numeric(), rt = numeric()),
  redundancy_filtering = TRUE
) {
  WGCNA::allowWGCNAThreads(parallel::detectCores())

  peaks <- peaks %>%
    dplyr::select(mz, rt, dplyr::starts_with("intensity")) %>%
    dplyr::mutate(peak = lexicographic_rank(mz, rt)) %>%
    dplyr::distinct(peak, .keep_all = TRUE) %>%
    tidyr::pack(intensity = dplyr::starts_with("intensity"))
  intensity_matrix <- t(dplyr::pull(peaks, intensity)) %>%
    magrittr::set_colnames(dplyr::pull(peaks, peak))
  correlation_matrix <- WGCNA::cor(intensity_matrix, use = "p", method = "p")

  clustering <- peak_clustering(
    intensity_matrix = intensity_matrix,
    correlation_matrix = correlation_matrix,
    correlation_threshold = correlation_threshold,
    min_cluster_size = min_cluster_size,
    network_type = network_type,
    deep_split = deep_split
  )

  peaks <- peaks %>%
    tibble::add_column(cluster = clustering) %>%
    dplyr::mutate(intensity = rowMeans(intensity, na.rm = TRUE)) %>%
    dplyr::mutate(mass_defect = mass_defect(mz, precision = 0.01))

  annotation <- simple_annotation(peaks, compounds, adducts, mass_tolerance) %>%
    dplyr::inner_join(peaks, by = c("mz", "rt"))

  annotation <- compute_scores(annotation, adduct_weights, time_tolerance)
  # annotation <- evaluate_pathways(annotation, pathways, score_threshold = 0.1)
  annotation <- compute_confidence_levels(annotation, expected_adducts) %>%
    boost_scores(boost, mass_tolerance, time_tolerance) %>%
    print_confidence_distribution()

  if (!redundancy_filtering)
    return(annotation)

  annotation <- annotation %>%
    remove_duplicates(score_threshold = 0) %>%
    print_confidence_distribution()
}
