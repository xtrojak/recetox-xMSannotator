.is_nonunique <- function(values) {
  frequency <- table(values)
  values %in% names(frequency)[frequency > 1]
}

.annotate_by_ovelapping_mz <- function(peak, meta, mz_tolerance_ppm) {
  threshold <- as.double(mz_tolerance_ppm) / 1000000

  indices <- tidyr::expand_grid(i=seq_len(nrow(peak)), j=seq_len(nrow(meta)))
  indices <- dplyr::filter(indices, abs(peak$mz[i] - meta$mz[j]) <= peak$mz[i] * threshold)

  colnames(peak) <- paste0("peak.", colnames(peak))
  colnames(meta) <- paste0("metabolite.", colnames(meta))
  annotation <- purrr::pmap_dfr(indices, function(i, j) cbind(peak[i,], meta[j,]))

  # call garbage collector after heavy memory load
  rm(indices)
  gc()

  annotation
}

simple_annotation <- function(data, metabolites, mz_tolerance_ppm = 10, adduct_selection = NULL) {
  stopifnot(all(c('mz', 'rt') %in% colnames(data)))
  stopifnot(all(c('adduct', 'formula', 'mz') %in% colnames(metabolites)))

  data <- dplyr::distinct(data, mz, rt, .keep_all = TRUE)
  metabolites <- dplyr::distinct(metabolites, .keep_all = TRUE)

  adduct_selection <- as.character(adduct_selection)
  if (length(adduct_selection)) {
    metabolites <- metabolites[metabolites$adduct %in% adduct_selection,]
  }

  annotation <- .annotate_by_ovelapping_mz(data, metabolites, mz_tolerance_ppm)
  #annotation <- dplyr::distinct(annotation, .keep_all = TRUE)

  if (!nrow(annotation)) {
    print("no matches found")
    return(annotation)
  }

  #annotation <- .remove_invalid_formulas(annotation)
  #annotation <- .remove_invalid_water_adducts(annotation)
  #annotation$multiple_match <- .is_nonunique(annotation$mz)

  return(as.data.frame(annotation))
}
