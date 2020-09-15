print_confidence_distribution <- function(df) {
  print("Stage 4 confidence level distribution for unique metabolites")
  print(df %>% filter(!duplicated(metabolite)) %>% count(confidence))

  print("Stage 4 confidence level distribution for unique formulas")
  print(df %>% filter(!duplicated(formula)) %>% count(confidence))

  invisible(df)
}

redundancy_filtering <- function(annotation, score_threshold) {
  is_redundant <- function (adduct_weight, score) {
    score <- ifelse(is.na(adduct_weight), score, 100 * score)
    score != max(score)
  }

  recompute_score <- function (redundant, score) {
    ((length(redundant) - sum(redundant)) / length(redundant)) * score
  }

  annotation %>%
    group_by(mz) %>%
    mutate(redundant = is_redundant(adduct_weight, score)) %>%
    group_by(metabolite) %>%
    mutate(score = recompute_score(redundant, score)) %>%
    ungroup() %>%
    filter(!redundant, score >= {{ score_threshold }}) %>%
    select(-redundant, -multiple_match) %>%
    print_confidence_distribution()
}