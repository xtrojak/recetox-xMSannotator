evaluate_pathways <- function(annotation, pathways, score_threshold) {
  # TODO: exluded_pathways <- "map01100"
  # TODO: exluded_metabolites <- c("HMDB29244", "HMDB29245", "HMDB29246")

  pathways <- distinct(pathways, pathway, metabolite)
  annotation <- left_join(annotation, pathways, by = "metabolite")

  # TODO: filter metabolites
  selection <- NULL

  max_cardinality <- annotation %>%
    select(metabolite, cluster) %>%
    add_count(cluster) %>%
    group_by(metabolite) %>%
    transmute(n = max(n)) %>%
    pull()

  mutate(
    annotation,
    score = ifelse(metabolite %in% {{ selection }}, score + {{ max_cardinality }}, score)
  )
}

#function (df, cluster, pathway) {
#  g <- df$pathway %in% pathway
#  f <- df$cluster %in% cluster
#
#  a <- n_distinct(df$metabolite[f & g])
#  b <- n_distinct(df$metabolite[f & !g])
#  c <- n_distinct(df$metabolite[!f & g])
#  d <- n_distinct(df$metabolite[!f & !g])
#
#  fisher.test(matrix(a, c, b, d, nrow = 2))$p.value
#}


#selection <- annotation %>%
#  filter(score > {{ score_threshold }} & !is.na(adduct_weight)) %>%
#  pull(metabolite)
#distinct(pathways, pathway, metabolite) %>%
#  mutate(
#    selected = metabolite %in% {{ selection }},
#    b = n_distinct(metabolite[selected]),
#    d = n_distinct(metabolite[!selected]),
#  ) %>%
#  group_by(pathway) %>%
#  summarise(
#    a = n_distinct(metabolite[selected]),
#    c = n_distinct(metabolite[!selected]),
#    b = b - a,
#    d = d - c,
#  ) %>% left_join()
#  filter(fisher.test(matrix(a, c, b, d, nrow = 2))$.p.value <= 0.05)