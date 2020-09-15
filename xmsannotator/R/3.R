# annotation
#   metabolite
#   adduct
#   score
#   peak_cluster
#
# pathways
#   pathway
#   metabolite


selection_independence_on_pathway_association <- function(df, score_threshold) {
  df %>%
    mutate(
      selection = replace_na(score, 0) >= {{ score_threshold }} & !is.na(adduct_weight),
      total_selection = n_distinct(metabolite[selection]),
      total_nselection = n_distinct(metabolite[!selection])) %>%
    group_by(pathway) %>%
    mutate(selection_independence_on_pathway_association = {
      a <- n_distinct(metabolite[selection])
      b <- total_selection - a
      c <- n_distinct(metabolite[!selection])
      d <- total_nselection - c
      fisher.test(matrix(a, b, c, d, nrow = 2))$p.value <= 0.05
    }) %>%
    ungroup()
}

evaluate_pathways <- function(annotation, pathways, score_threshold) {
  exluded_pathways <- "map01100"
  exluded_metabolites <- c("HMDB29244", "HMDB29245", "HMDB29246")

  annotation %>% select(metabolite, score, adduct_weight)

  pathways %>%
    distinct(pathway, metabolite) %>%
    left_join(, by=metabolite)%>%
    selection_independence_on_pathway_association(score_threshold) %>%
    filter(selection_independence_on_pathway_association & selection)

  filter(pathways, {})
  mutate(annotation, score = ifelse())
}
