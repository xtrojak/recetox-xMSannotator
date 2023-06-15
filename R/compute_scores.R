#' Compute a confidence score for each annotation.
#'
#' @param annotation A table with annotated peaks.
#' @param adduct_weights A weight-by-adduct table.
#'
#' @import dplyr
#' @importFrom rlang .data
compute_scores <- function(annotation,
                           adduct_weights) {
  annotation <- left_join(annotation, adduct_weights, by = "adduct")
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, score = 10^max(-Inf, weight, na.rm = TRUE))
  annotation <- ungroup(annotation)
  annotation <- compute_multiple_matches(annotation)
  annotation
}
