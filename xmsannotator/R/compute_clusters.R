peak_clustering <- function(intensity_matrix, correlation_matrix, correlation_threshold, deep_split, min_cluster_size, network_type) {
  inverse_correlation_matrix <- 1 - correlation_matrix
  power_estimate <- WGCNA::pickSoftThreshold(intensity_matrix)$powerEstimate

  if (is.na(power_estimate)) {
    warning("Unable to estimate soft threshold. Using fallback value.")
    power_estimate <- 6
  }

  # NOTE: When WGCNA is loaded incorectly then WGCNA::blockwiseModules may fail
  # due to conflicting functions WGCNA::cor and stats::cor.
  modules <- WGCNA::blockwiseModules(
    datExpr = intensity_matrix,
    checkMissingData = FALSE,
    blocks = dynamicTreeCut::cutreeDynamic(
      dendro = flashClust::flashClust(as.dist(inverse_correlation_matrix)),
      distM = inverse_correlation_matrix,
      minClusterSize = min_cluster_size,
      pamRespectsDendro = FALSE
    ),
    blockSizePenaltyPower = 100,
    power = power_estimate,
    deepSplit = deep_split,
    detectCutHeight = NULL,
    minModuleSize = min_cluster_size,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 1 - correlation_threshold
  )

  if (!modules$MEsOK) {
    warning("Unable to perform peak WGCNA clustering. Using fallback method.")

    distance_matrix <- WGCNA::TOMdist(
      WGCNA::adjacency(
        datExpr = intensity_matrix,
        type = network_type,
        power = power_estimate,
        corOptions = list(use = "p", method = "pearson")
      )
    )

    modules <- WGCNA::mergeCloseModules(
      exprData = intensity_matrix,
      colors = dynamicTreeCut::cutreeDynamic(
        dendro = flashClust::flashClust(as.dist(distance_matrix)),
        distM = distance_matrix,
        deepSplit = deep_split,
        minClusterSize = min_cluster_size,
        pamRespectsDendro = FALSE
      ),
      cutHeight = 1 - correlation_threshold,
      trapErrors = TRUE
    )
  }

  as.integer(as.factor(modules$colors))
}