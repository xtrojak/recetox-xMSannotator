.find_overlpapping_mzs <- function(data_a, data_b, mz_threshold) {
  ppmb <- mz_threshold * data_a$mz / 1000000

  result <- furrr::future_imap_dfr(data_a$mz, function(i, mz_a, mz_b, threshold) {
    indices <- which(abs(mz_b - mz_a) <= threshold)
    data.frame(index_a = rep(i, length(indices)), index_b = indices)
  }, mz_b = data_b$mz, threshold = ppmb)

  return(result)
}

.remove_water_adducts <- function(annotation, water_adducts = c("M+H-H2O", "M+H-2H2O", "M-H2O-H")) {
  is_water_adduct <- annotation$adduct %in% water_adducts

  if (any(is_water_adduct)) {
    is_valid_formula <- sapply(seq_along(is_water_adduct), function(i) {
      !is_water_adduct[i] || check_element(annotation$formula[i], elementname = "O") > 0
    })
    annotation <- annotation[is_valid_formula]
  }

  return(annotation)
}

.remove_invalid_formulas <- function (annotation, nops_check) {
  is_valid_formlula <- furrr::future_map_lgl(annotation$formula, .check_golden_rules, nops_check=nops_check)
  annotation <- annotation[is_valid_formlula,]

  return(annotation)
}

.check_golden_rules <- function(formula, nops_check) {
  formula <- as.character(formula)

  numnitrogens <- check_element(formula, "N")
  numcarbons <- check_element(formula, "C")
  numoxygens <- check_element(formula, "O")
  numhydrogens <- check_element(formula, "H")
  numphos <- check_element(formula, "P")
  numsulphur <- check_element(formula, "S")

  if (numcarbons < 1) {
    return(FALSE)
  }

  nitrogen_to_carbon_ratio <- numnitrogens / numcarbons
  oxygen_to_carbon_ratio <- numoxygens / numcarbons
  phosphorus_to_carbon_ratio <- numphos / numcarbons
  sulphur_to_carbon_ratio <- numsulphur / numcarbons
  hydrogens_to_carbon_ratio <- numhydrogens / numcarbons

  if (hydrogens_to_carbon_ratio < 0.1 |
    hydrogens_to_carbon_ratio > 6 |
    nitrogen_to_carbon_ratio > 4 |
    oxygen_to_carbon_ratio > 3 |
    phosphorus_to_carbon_ratio > 2 |
    sulphur_to_carbon_ratio > 3)
    return(FALSE)

  if (nops_check) {
    #NOPS>1
    if ((numnitrogens > 1 & numoxygens > 1 & numphos > 1 & numsulphur > 1) &
      (numnitrogens > 10 | numoxygens > 20 | numphos > 4 | numsulphur > 3))
      return(FALSE)

    #NOP>3
    if ((numnitrogens > 3 & numoxygens > 3 & numphos > 3) & (numnitrogens > 11 | numoxygens > 22 | numphos > 6))
      return(FALSE)

    #OPS>1
    if ((numoxygens > 1 & numphos > 1 & numsulphur > 1) & (numoxygens > 14 | numphos > 3 | numsulphur > 3))
      return(FALSE)

    #PSN>1
    if ((numnitrogens > 1 & numphos > 1 & numsulphur > 1) & (numnitrogens > 4 | numphos > 3 | numsulphur > 3))
      return(FALSE)

    #NOS>6
    if ((numnitrogens > 6 & numoxygens > 6 & numsulphur > 6) &
      (numnitrogens > 19 | numoxygens > 14 | numsulphur > 8))
      return(FALSE)
  }

  return(TRUE)
}

.is_nonunique_value <- function(values) {
  frequency <- table(values)
  values %in% names(frequency)[frequency > 1]
}

simple_annotation <- function(data, metabolite_db, max_mz_diff = 10, adduct_selection = NULL, workers = future::availableCores()) {
  future::plan(future::multiprocess(workers = workers))

  data <- unique(as.data.frame(data), by = c('mz', 'rt'))
  metabolite_db <- as.data.frame(metabolite_db)
  adduct_selection <- as.character(adduct_selection)

  metabolite_db <- unique(metabolite_db, by = c('formula', 'chemical_id'))
  metabolite_db <- metabolite_db[order(metabolite_db$chemical_id), ]

  if (length(adduct_selection)) {
    metabolite_db <- metabolite_db[metabolite_db$adduct %in% adduct_selection]
  }

  overlapping <- .find_overlpapping_mzs(data_a = data, data_b = metabolite_db, mz_threshold = max_mz_diff)
  annotation <- cbind(metabolite_db[overlapping$index_b,], data[overlapping$index_a,])
  annotation <- unique(annotation)

  if (!nrow(annotation)) {
    print("no matches found")
    return(annotation)
  }

  annotation <- .remove_invalid_formulas(annotation)
  annotation <- .remove_water_adducts(annotation)
  annotation$multiple_match <- .is_nonunique_value(annotation$mz)

  return(annotation)
}
