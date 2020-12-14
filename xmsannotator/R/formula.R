check_element <- function(formula, element) {
  stopifnot(length(element) == 1)

  pattern <- paste0(element, '(?![a-z]+)([0-9]*)')
  matches <- stringr::str_match_all(formula, pattern)

  sapply(matches, function(match) {
    counts <- as.integer(match[, 2])
    sum(counts, na.rm = TRUE) + sum(is.na(counts))
  })
}

check_golden_rules <- function(formula) {
  c <- check_element(formula, "C")
  h <- check_element(formula, "H")
  f <- check_element(formula, "F")
  n <- check_element(formula, "N")
  o <- check_element(formula, "O")
  p <- check_element(formula, "P")
  s <- check_element(formula, "S")

  ratio_check <- (0.1 <= h / c) &
    (h / c <= 6) &
    (f / c <= 6) &
    (n / c <= 4) &
    (o / c <= 3) &
    (p / c <= 2) &
    (s / c <= 3)
  ratio_check <- ifelse(c > 0, ratio_check, FALSE)

  nops <- !(n > 1 & o > 1 & p > 1 & s > 1) | (n < 10 & o < 20 & p < 4 & s < 3)
  nop <- !(n > 3 & o > 3 & p > 3) | (n < 11 & o < 22 & p < 6)
  ops <- !(o > 1 & p > 1 & s > 1) | (o < 14 & p < 3 & s < 3)
  psn <- !(n > 1 & p > 1 & s > 1) | (n < 4 & p < 3 & s < 3)
  nos <- !(n > 6 & o > 6 & s > 6) | (n < 19 & o < 14 & s < 8)

  return(ratio_check & nops & nop & ops & psn & nos)
}
