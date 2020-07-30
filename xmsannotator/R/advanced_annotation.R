.get_peak_blocks_modulesvhclust <- function(data, time_step, max_rt_diff, deep_split, min_cluster_size, cutheight, cormethod, networktype = "unsigned", workers, mycl_metabs = NA) {
  data_m <- t(data[, -(1:2)])
  colnames(data_m) <- paste(data$mz, data$rt, sep = "_")

  power_val <- tryCatch(expr = WGCNA::pickSoftThreshold(data = data_m, dataIsExpr = TRUE, powerVector = c(1:10, seq(from = 12, to = 20, by = 2)), verbose = 0)$powerEstimate,
                        error = function(...) return(6))

  mod_list <- try({
                    network <- WGCNA::blockwiseModules(datExpr = data_m, checkMissingData = FALSE, blocks = mycl_metabs, maxBlockSize = 5000, blockSizePenaltyPower = 100, randomSeed = 12345, loadTOM = FALSE, corType = "pearson", maxPOutliers = 1, quickCor = 0, pearsonFallback = "individual", cosineCorrelation = FALSE, power = power_val, networkType = "unsigned", TOMType = "signed", TOMDenom = "min", getTOMs = NULL, saveTOMs = FALSE, saveTOMFileBase = "blockwiseTOM", deepSplit = deep_split, detectCutHeight = NULL, minModuleSize = min_cluster_size, maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, minSplitHeight = NULL, minAbsSplitHeight = NULL, useBranchEigennodeDissim = FALSE, minBranchEigennodeDissim = mergeCutHeight, stabilityLabels = NULL, minStabilityDissim = NULL, pamStage = TRUE, pamRespectsDendro = FALSE, reassignThreshold = 1e-06, minCoreKME = 0.5, minCoreKMESize = min_cluster_size / 3, minKMEtoStay = 0.3, mergeCutHeight = cutheight, impute = TRUE, trapErrors = FALSE, numericLabels = FALSE, nThreads = workers, verbose = 0, indent = 0)
                    as.numeric(network$colors)
                  }, silent = TRUE)

  if (inherits(mod_list, 'try-error')) {
    adjacency_matrix <- WGCNA::adjacency(datExpr = data_m, type = networktype, power = power_val, corOptions = list(use = 'p', method = cormethod))
    if (any(duplicated_names <- duplicated(rownames(adjacency_matrix)))) {
      adjacency_matrix <- adjacency_matrix[!duplicated_names, !duplicated_names]
    }

    dissTOMCormat <- WGCNA::TOMdist(adjacency_matrix)
    hierTOMCormat <- WGCNA::flashClust(as.dist(dissTOMCormat), method = "complete")

    mod_list <- dynamicTreeCut::cutreeDynamic(hierTOMCormat, distM = dissTOMCormat, deepSplit = deep_split, minClusterSize = min_cluster_size, pamRespectsDendro = FALSE, pamStage = TRUE)
    if (length(unique(mod_list)) > 1) {
      try(mod_list <- as.numeric(WGCNA::mergeCloseModules(data_m, colors = colorhdataOne2, cutHeight = cutheight)$colors), silent = TRUE)
    }
  }

  t1 <- table(mod_list)
  mod_names <- names(t1)
  mod_names <- as.numeric(mod_names)

  diffmatB <- purrr::map_dfr(mod_names, function(group_number) {
    subdata <- data[mod_list == group_number,]
    subdata <- subdata[order(subdata$rt),]

    groupB <- group_by_rt_histv1(subdata, time_step, max_diff_rt = max_rt_diff, groupnum = group_number)
    rownames(groupB) <- NULL

    return(groupB)
  })

  return(diffmatB)
}

.clusterring <- function(data, correlation_method, correlation_threshold, min_cluster_size, deep_split, network_type, workers) {
  stopifnot(c('mz', 'rt') %in% colnames(data))

  correlation <- WGCNA::cor(x = t(data[, -(1:2)]), nThreads = workers, method = correlation_method, use = 'p')
  dendrogram <- flashClust::flashClust(as.dist(1 - correlation), method = 'complete')
  mycl_metabs <- dynamicTreeCut::cutreeDynamic(dendro = dendrogram, distM = 1 - correlation, deepSplit = 1, minClusterSize = min_cluster_size, pamStage = TRUE, pamRespectsDendro = FALSE, verbose = 0)

  clustering <- .get_peak_blocks_modulesvhclust(dataA = data, time_step = 1, max.rt.diff = max_diff_rt, cor.thresh = NA, deep_split = deep_split, minclustsize = min_cluster_size, cutheight = 1 - correlation_threshold, cormethod = correlation_method, networktype = network_type, workers = workers, mycl_metabs = mycl_metabs)
  clustering <- dplyr::arrange(clustering, mz, rt)

  return(clustering)
}

.step5 <- function(result, adduct_weights, score_threshold = 0) {
  result <- dplyr::arrange(result, dplyr::desc(confidence), dplyr::desc(chemical_id), dplyr::desc(score), dplyr::desc(adduct))

  mz_frequency <- table(result$mz)
  nonunique_mzs <- names(mz_frequency)[mz_frequency > 1]

  bad_indices <- purrr::map_int(nonunique_mzs, function(value) {
    indices <- which(result$mz %in% value)
    subset <- result[indices,]

    good_adduct <- which(subset$adduct %in% adduct_weights$adduct)
    subset$score[good_adduct] <- subset$score[good_adduct] * 100

    indices <- indices[which(subset$score != max(subset$score, na.rm = TRUE))]

    for (chemical_id in subset$chemical_id[indices]) {
      ind <- which(result$chemical_id %in% chemical_id)
      result$score[ind] <- (length(ind) - 1) * result$score[ind[1]] / length(ind)
    }

    indices
  })

  result <- result[-bad_indices]
  result <- result[result$score >= score_threshold,]
  result$multiple_match <- .is_nonunique_value(result$mz)

  print("Stage 5 confidence level distribution for unique chemical/metabolite IDs")
  print(table(unique(result, by = 'chemical_id')$confidence))

  print("Stage 5 confidence level distribution for unique chemical/metabolite formulas")
  print(table(unique(result, by = 'formula')$confidence))

  return(result)
}

advanced_annotation <- function(data, metabolites, adduct_weights, max_mz_diff, max_rt_diff, correlation_method, correlation_threshold, min_cluster_size, deep_split, network_type, boost_metabolites, expected_adducts, min_isp, max_isp, strict_boosting, redundancy_filtering, workers = future::availableCores()) {
  data <- tibble::as_tibble(data)
  metabolites <- tibble::as_tibble(metabolites)
  adduct_weights <- tibble::as_tibble(adduct_weights)
  max_mz_diff <- as.double(max_mz_diff)
  max_rt_diff <- as.double(max_rt_diff)
  correlation_method <- as.character(correlation_method)
  correlation_threshold <- as.double(correlation_threshold)
  min_cluster_size <- as.integer(min_cluster_size)
  deep_split <- as.integer(deep_split)
  network_type <- as.character(network_type)
  boost_metabolites <- as.data.frame(boost_metabolites)
  expected_adducts <- as.character(expected_adducts)
  min_isp <- as.integer(min_isp)
  max_isp <- as.integer(max_isp)
  strict_boosting <- as.logical(strict_boosting)
  redundancy_filtering <- as.logical(redundancy_filtering)

  stopifnot(all(c('mz', 'rt') %in% colnames(data)))
  stopifnot(all(c("mz", "name", "formula", "monoisotopic_mass", "adduct", "adduct_mass") %in% colnames(metabolites)))

  data <- dplyr::distinct(data, mz, rt, .keep_all = TRUE)

  if (is.null(adduct_weights)) {
    adduct_weights <- tibble::tibble(Adduct = c('M+H', 'M-H'), Weight = c(1, 1))
  }

  data(adduct_table)
  outloc <- getwd()

  xMSannotator::multilevelannotation(
    dataA = data, outloc = outloc, max.mz.diff = max_mz_diff, max.rt.diff = max_rt_diff,
    cormethod = correlation_method, corthresh = correlation_threshold, db_name = 'Custom',
    adduct_table = adduct_table, adduct_weights = adduct_weights, num_nodes = workers,
    customIDs = NA, deepsplit = deep_split, networktype = network_type,
    minclustsize = min_cluster_size, redundancy_check = redundancy_filtering,
    filter.by = expected_adducts, min_ions_perchem = min_isp, boostIDs = boost_metabolites,
    max_isp = max_isp, customDB = metabolites, pathwaycheckmode = if (strict_boosting) 'pm' else 'p'
  )

  # future::plan(future::multiprocess(workers = workers))
  # stage_01 <- .clusterring(data = data, correlation_method = correlation_method, correlation_threshold = correlation_threshold, min_cluster_size = min_cluster_size, deep_split = deep_split, network_type = network_type, workers = workers)
  # stage_02 <- simple_annotation(data = data[, c('mz', 'rt')], metabolite_db = metabolites_db, adduct_selection = NULL, max_mz_diff = max_mz_diff, workers = workers)
  # return('stage_01' = stage_01, 'stage_02' = stage_02)
}
