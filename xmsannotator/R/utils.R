lexicographic_rank <- function (...) {
  .o <- order(...)
  .x <- cbind(...)[.o, ]
  cumsum(!duplicated(.x))[order(.o)]
}

save_parquet <- function(data, file) {
  invisible(arrow::write_parquet(data, file))
}

load_parquet <- function(file, required_columns, optional_columns = NULL) {
  df <- tryCatch(
    arrow::read_parquet(file, col_select = dplyr::any_of(union(required_columns, optional_columns))),
    error = function(e) stop("The file ", toString(file), " is not a valid Parquet file")
  )

  if (!all(required_columns %in% colnames(df))) {
    stop("The file ", toString(file), " does not contain any of these required columns: ", toString(required_columns))
  }

  return(df)
}

load_csv <- function (file, required_columns, optional_columns = NULL) {
  df <- readr::read_csv(file)
  df <- dplyr::select(df, dplyr::any_of(union(required_columns, optional_columns)))

  if (!all(required_columns %in% colnames(df))) {
    stop("The file ", toString(file), " does not contain any of these required columns: ", toString(required_columns))
  }

  return(df)
}

load_peak_table_parquet <- function(file) {
  load_parquet(file, required_columns = c("peak", "mz", "rt"))
}

load_adduct_table_parquet <- function(file) {
  load_parquet(file, required_columns = c("adduct", "charge", "mass", "n_molecules"))
}

load_compound_table_parquet <- function(file) {
  compound_table <- load_parquet(
    file,
    required_columns = c("monoisotopic_mass", "molecular_formula"),
    optional_columns = c("compound", "recetox_cid")
  )

  if (!any(c("compound", "recetox_cid") %in% colnames(compound_table))) {
    stop("The supplied compound table needs to contain either compound or recetox_cid column")
  }

  if ("recetox_cid" %in% colnames(compound_table)) {
    compound_table <- dplyr::mutate(compound_table, compound = recetox_cid)
  }

  return(compound_table)
}

load_peak_table_hdf <- function(file) {
  peak_table <- rhdf5::h5read(file, "peaks")
  peak_table <- dplyr::select(
    peak_table,
    dplyr::any_of("peak"),
    dplyr::all_of(c("mz", "rt")),
    dplyr::starts_with("intensity")
  )

  if (!"peak" %in% colnames(peak_table)) {
    peak_table <- dplyr::mutate(peak_table, peak = lexicographic_rank(mz, rt))
  }

  return(peak_table)
}

load_expected_adducts_csv <- function (file) {
  load_csv(file, required_columns = "adduct")
}

load_boost_compounds_csv <- function (file) {
  load_csv(file, required_columns = "compound", optional_columns = c("mz", "rt"))
}
