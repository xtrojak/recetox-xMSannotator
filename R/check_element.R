#' @export
check_element <- function(curformula, elementname) {
  curformula <- as.character(curformula)
  formula_split <- strsplit(curformula, split = "")
  
  element_length <- length(elementname)
  
  match <- gregexpr(curformula, pattern = paste(elementname, "[0-9]*", sep = ""))
  
  numelement <- 0
  
  if (length(match) > 0) {
    match_length <- attr(match[[1]], "match.length")[1] - element_length 
    if (element_length > 1 && (match_length + element_length > 2)) {
      match_length <- match_length + 1
    }
    
    if (match_length > 0) {
      numelement <- paste(formula_split[[1]][(match[[1]][1] + element_length):(match[[1]][1] + match_length)], collapse = "")
      numelement <- as.numeric(numelement)
    } else if (match_length == 0) {
      numelement <- 1
    }
  }
  return(numelement)
}
