## recetox-aplcms to xMSannotator

### R Sample Code

``` R
#' @param file Path to parquet file to be loaded into memory
#'
#' @export
load_parquet <- function(file) {
    message <- paste("The file", toString(file), "seams not to be a valid Parquet file.")
    rlang::with_handlers(
        arrow::read_parquet(file),
        error = ~ rlang::abort(message, parent = .)
    )
}


#' @param peak_table Peak table from recetox.aplcms package.
#'
#' @export
#' @import dplyr
rcx_aplcms_to_xmsannotator <- function(peak_table) {
    col_base <- c("feature", "mz", "rt")
    output_table <- peak_table %>% distinct(across(any_of(col_base)))

    for (level in levels(peak_table$sample)) {
        subdata <- peak_table %>%
            filter(sample == level) %>%
            select(any_of(c(col_base, "sample_intensity"))) %>%
            rename(!!level := "sample_intensity")
        output_table <- inner_join(output_table, subdata, by = col_base)
    }

    output_table <- output_table  %>% select(-feature) %>% rename(time = rt)
    return(output_table)
}

```

## recetox-aplcms to recetox-xMSannotator
This is the same as to the recetox-xMSannotator, but the `feature` column has to be renamed to peak, and `rt` can stay the same.

### R Sample Code

```R
#' @param peak_table Peak table from recetox.aplcms package.
#'
#' @export
#' @import dplyr
rcx_aplcms_to_rcx_xmsannotator <- function(peak_table) {
    col_base <- c("feature", "mz", "rt")
    output_table <- peak_table %>% distinct(across(any_of(col_base)))

    for (level in levels(peak_table$sample)) {
        subdata <- peak_table %>%
            filter(sample == level) %>%
            select(any_of(c(col_base, "sample_intensity"))) %>%
            rename(!!level := "sample_intensity")
        output_table <- inner_join(output_table, subdata, by = col_base)
    }

    output_table <- output_table %>% rename(peak = feature)
    return(output_table)
}

```

## HMDB: Rda to Parquet
These are simple instructions how to convert `RDA` file containing the HMDB annotation database (used in xMSannotator) to `parquet` format.

The important steps are defining `filter_keys` and `target_keys`, which are vectors containing keywords to be selected and subsequently renamed, in the respective order.

### R Sample Code

```R
save_parquet <- function(data, filename) {
  invisible(arrow::write_parquet(data, filename))
}

load_rda <- function(filename) {
  get(load(filename))
}

convert_rda_to_parquet <- function(rda, filter_keys, target_keys) {
  target_data <- rda[filter_keys]
  names(target_data) <- target_keys
  return(target_data)
}

rda <- load_rda("filename.rda")
filter_keys <- c("Name", "InChiKey", "Formula", "MonoisotopicMass")
target_keys <- c("iupac_name", "iupac_inchikey", "molecular_formula", "monoisotopic_mass")
result <- convert_rda_to_parquet(rda, filter_keys, target_keys)
save_parquet(result, "filename.parquet")

```