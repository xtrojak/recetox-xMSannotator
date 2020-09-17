validate_peak_table <- function(peaks) {
  columns <- c("mz", "rt")

  if (!all(columns %in% colnames(peaks)))
    stop("Provided peak table must contain these columns: ", toString(columns))

  distinct(peaks, mz, rt, .keep_all = TRUE)
}


validate_adduct_table <- function(adducts) {
  columns <- c("adduct", "charge", "mass", "molecules")

  if (!all(columns %in% colnames(adducts)))
    stop("Provided adduct table must contain these columns: ", toString(columns))
  if (anyDuplicated(adducts$adduct))
    stop("Provided adduct table contains duplicated adducts!")

  select(adducts, all_of(columns))
}


validate_metabolite_table <- function(metabolites) {
  columns <- c("metabolite", "monoisotopic_mass", "formula")

  if (!all(columns %in% colnames(metabolites)))
    stop("Provided metabolite table must contain these columns: ", toString(columns))
  if (anyDuplicated(metabolites$metabolite))
    stop("Provided metabolite table contains duplicated metabolites!")

  select(metabolites, all_of(columns))
}


compile_metabolites <- function(adducts, metabolites) {
  expand_grid(metabolites, adducts) %>%
    mutate(expected_mass = (molecules * monoisotopic_mass + mass) / charge) %>%
    select(all_of(union(colnames(metabolites), c("adduct", "expected_mass"))))
}


merge_by_mz <- function(a, b, mz_tolerance) {
  a_indices <- seq_len(nrow(a))
  b_indices <- seq_len(nrow(b))

  indices <- expand_grid(i = a_indices, j = b_indices) %>%
    filter(abs(a$mz[i] - b$expected_mass[j]) <= a$mz[i] * mz_tolerance)
  gc() # call garbage collector after heavy memory load

  a <- slice(a, indices$i)
  b <- slice(b, indices$j)

  bind_cols(a, b)
}


is_nonuinque_mz <- function(mz) {
  # FIXME: use grouping with parametric tolerance (use rounding if nothing else)
  as_tibble_col(mz) %>% add_count(value) %>% transmute(n > 1) %>% pull()
}


simple_annotation <- function(peaks, metabolites, adducts, mz_tolerance_ppm = 10) {
  peaks <- validate_peak_table(peaks)
  adducts <- validate_adduct_table(adducts)
  metabolites <- validate_metabolite_table(metabolites)
  metabolites <- compile_metabolites(adducts, metabolites)

  merge_by_mz(peaks, metabolites, mz_tolerance_ppm / 10^6) %>%
    filter(check_golden_rules(formula)) %>%
    filter(is_valid_adduct(adduct, formula)) %>%
    mutate(multiple_match = is_nonuinque_mz(mz)) %>%
    as.data.frame()
}
