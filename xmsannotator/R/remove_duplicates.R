remove_duplicates <- function(annotation, score_threshold) {
  is_redundant <- function(weight, score) {
    score <- score * 100^is.na(weight)
    score != max(score)
  }

  recompute_score <- function(redundant, score) {
    ((length(redundant) - sum(redundant)) / length(redundant)) * score
  }

  annotation %>%
    dplyr::group_by(mz) %>%
    dplyr::mutate(redundant = is_redundant(adduct_weight, score)) %>%
    dplyr::group_by(compound) %>%
    dplyr::mutate(score = recompute_score(redundant, score)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!redundant, score >= { { score_threshold } }) %>%
    dplyr::select(-redundant, -multiple_match)
}
