# peak_intensity_table
#   mz
#   rt
#   intensity_<file>
# -----
#   module
#   rt_cluster
#   is_cluster
#   mz_defect
#   multiple_match
#   score

# metabolites
#   mz
#   chemical_id
#   name
#   formula
#   monoisotopic_mass
#   adduct
#   adduct_mass

# adduct_weights
#   adduct
#   adduct_weight


nest_clusters <- function(...) {
  factor(paste(..., sep = "_"))
}

cluster_by_density <- function(x, tolerance) {
  # FIXME: Though this implementation is same as the original implementation of
  # group_by_histv1, I dont thing that this clustering algorithm is correct or
  # even suitable for the clustering data (mainly retention time data).
  # For the future, consider something as dbscan which is also a density based
  # clustering algorithm.
  # Note, that the original condition `... && abs(h$mid - minimum) > tolerance`
  # smells bad.

  bandwith <- try(min(bw.nrd(x), tolerance))
  if (inherits(bandwith, "try-error"))
    bandwith <- tolerance

  bins <- seq(min(x) - tolerance, max(x) + tolerance, bandwith)
  histogram <- hist(x, breaks = bins, plot = FALSE)

  counter <- 1
  labels <- integer(length(histogram$density))
  edges <- histogram$density == 0 & abs(histogram$mids - min(x)) > tolerance

  for (i in seq_along(labels)) {
    if (edges[i])
      counter <- counter + 1
    labels[i] <- if (histogram$density[i] > 0) counter else 0
  }

  selection <- histogram$density > 0
  labels <- labels[selection]
  mids <- histogram$mids[selection]
  mids[1] <- min(x)

  sapply(x, function(value) {
    labels[which.min(abs(mids - value))]
  })
}

peak_intensity_clustering <- function(peaks, peak_intensity_matrix, peak_correlation_matrix, correlation_threshold, deep_split, min_cluster_size, network_type) {
  power <- try(WGCNA::pickSoftThreshold(peak_intensity_matrix)$powerEstimate)
  if (is.na(power) || inherits(power, "try-error")) {
    power <- 6
  }

  dendrogram <- flashClust::flashClust(as.dist(1 - peak_correlation_matrix))
  dendrogram <- dynamicTreeCut::cutreeDynamic(
    dendro = dendrogram,
    distM = 1 - peak_correlation_matrix,
    minClusterSize = min_cluster_size,
    pamRespectsDendro = FALSE
  )

  modules <- try(WGCNA::blockwiseModules(
      datExpr = peak_intensity_matrix,
      checkMissingData = FALSE,
      blocks = dendrogram,
      blockSizePenaltyPower = 100,
      power = power,
      deepSplit = deep_split,
      detectCutHeight = NULL,
      minModuleSize = min_cluster_size,
      pamRespectsDendro = FALSE,
      mergeCutHeight = 1 - correlation_threshold
    )$colors)

  if (inherits(modules, "try-error")) {
    adjacency_matrix <- WGCNA::adjacency(
      datExpr = peak_intensity_matrix,
      type = network_type,
      power = power,
      corOptions = list(use = "p", method = "pearson")
    )
    distance_matrix <- WGCNA::TOMdist(adjacency_matrix)

    dendrogram <- flashClust::flashClust(
      as.dist(distance_matrix)
    )
    dendrogram <- dynamicTreeCut::cutreeDynamic(
      dendro = dendrogram,
      distM = distance_matrix,
      deepSplit = deep_split,
      minClusterSize = min_cluster_size,
      pamRespectsDendro = FALSE
    )
    modules <- WGCNA::mergeCloseModules(
      exprData = peak_intensity_matrix,
      colors = dendrogram,
      cutHeight = 1 - correlation_threshold,
      trapErrors = TRUE
    )$colors
  }

  tibble::add_column(peaks, peak_cluster = modules)
}

retention_time_clustering <- function(df, rt_tolerance) {
  as_tibble(df) %>%
    group_by(peak_cluster) %>%
    mutate(
      rt_cluster = cluster_by_density(mz, {{ rt_tolerance }}),
      rt_cluster = nest_clusters(peak_cluster, rt_cluster)
    ) %>%
    ungroup(peak_cluster)
}

mass_defect_clustering <- function(df, precision) {
  mutate(df,
    mass_defect = cut(mz %% 1, seq(0, 1, {{ precision }}), labels = FALSE),
    mass_defect_cluster = nest_clusters(rt_cluster, mass_defect),
  )
}

compute_delta_mz <- function(df) {
  # NOTE: I changed the computation of delta_ppm becaouse I think the previos
  # version behaved incorectly (see lines 130-140 in multilevelannotationstep4).
  mutate(df, delta_ppm = round(10^6 * abs(expected_mass - mz) / expected_mass, 2))
}

#compute_scores <- function(df, corelation_threshold, correlation_matrix, adduct_weights, expected_adducts, max_isp, rt_tolerance) {
#  data.frame(
#    mz = mz,
#    time = rt,
#    MatchCategory = ifelse(multiple_match, "Multiple", "Unique"),
#    theoretical.mz =,
#    chemical_ID = chemical_id,
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
#    group_by(chemical_id) %>%
#    mutate(score = 0) %>%
#    ungroup()
#}

compute_scores <- function(annotation, adduct_weights) {
  annotation <- as_tibble(annotation) %>%
    filter(is_valid_adduct(adduct, formula)) %>%
    left_join(adduct_weights, by = "adduct")
  mutate(annotation, score = 0)
}

advanced_annotation <- function(peaks, metabolites, adducts, mz_tolerance_ppm = 10, rt_tolerance = 10, correlation_threshold = 0.7, deep_split = 2, min_cluster_size = 10, network_type = "unsigned", adduct_weights, boost_metabolites, expected_adducts, workers = parallel::detectCores()) {
  WGCNA::allowWGCNAThreads(workers)

  peaks <- distinct(peaks, mz, rt, .keep_all = TRUE)
  metabolites <- distinct(metabolites)

  peak_intensity_matrix <- t(select(peaks, starts_with("intensity")))
  peak_correlation_matrix <- WGCNA::cor(peak_intensity_matrix, use = "p", method = "p")

  peaks %>%
    peak_intensity_clustering(
      peak_intensity_matrix = peak_intensity_matrix,
      peak_correlation_matrix = peak_correlation_matrix,
      correlation_threshold = correlation_threshold,
      min_cluster_size = min_cluster_size,
      network_type = network_type,
      deep_split = deep_split
    ) %>%
    retention_time_clustering(rt_tolerance = rt_tolerance) %>%
    mass_defect_clustering(precision = 0.01) %>%
    xmsannotator::simple_annotation(
      adducts = adducts,
      metabolites = metabolites,
      mz_tolerance_ppm = mz_tolerance_ppm
    ) %>%
    compute_scores(adduct_weights = adduct_weights) %>%
    #compute_scores(
    #  correlation_threshold = correlation_threshold,
    #  correlation_matrix = correlation_matrix,
    #  adduct_weights = adduct_weights,
    #  expected_adducts = as.character(expected_adducts),
    #  max_isp = as.integer(max_isp),
    #  rt_tolerance = rt_tolerance,
    #)
    evaluate_pathways(pathways = , score_threshold = 0.1) %>%
    step4(
      expected_adducts = expected_adducts,
      boost_metabolites = boost_metabolites,
      mz_tolerance = 10^6 * mz_tolerance_ppm,
      rt_tolerance = rt_tolerance
    ) %>%
    redundancy_filtering(score_threshold = 0)
}
