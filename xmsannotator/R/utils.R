validate_table <- function(table, primary_key = NULL, expected_columns = NULL, optional_columns = NULL) {
  expected_columns <- union(primary_key, expected_columns)

  if (!all(expected_columns %in% colnames(table)))
    stop("The table must containt these columns: ", toString(expected_columns))

  table <- table %>%
    dplyr::select(all_of(expected_columns), any_of(optional_columns)) %>%
    dplyr::distinct()

  if (!is.null(primary_key) && anyDuplicated(dplyr::select(table, primary_key)))
    stop("Duplicate values in primary key detected: ", toString(primary_key))
  table
}

load_hdf <- function (path, key) {
  if (!is.character(path) || !file.exists(path) || !hdf5r::is_hdf5(path))
    stop("File \"", path, "\" is not a valid HDF5 file.")

  file <- hdf5r::H5File$new(filename = path, mode = "r")
  on.exit(file$close_all())

  if (!key %in% names(file))
    stop("File \"", path, "\" does not contains required key \"", key, "\".")

  tibble::as_tibble(file[[key]]$read())
}

save_hdf <- function(path, key, data) {
  file <- hdf5r::H5File$new(filename = path, mode = "w")
  on.exit(file$close_all())

  file[[as.character(key)]] <- data
  invisible(data)
}
