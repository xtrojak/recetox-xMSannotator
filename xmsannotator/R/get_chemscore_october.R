
#' @import tidyr
#' @import dplyr
#' @import plyr
#' @importFrom magrittr %>%
#' @export
get_chemscore <- function(...,
                          annotation,
                          corthresh,
                          global_cor,
                          max_diff_rt = 10,
                          adduct_weights,
                          filter.by = c("M+H"),
                          mass_defect_window = 0.01,
                          MplusH.abundance.ratio.check = TRUE,
                          outlocorig) {

  query <- tibble(...)
  setwd(outlocorig)

  outloc1 <- paste(outlocorig, "/stage2/", sep = "")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1)

  curmchemdata <- dplyr::filter(
    annotation,
    chemical_ID == query$chemical_ID &
    abs(time - query$time) <= 10
  )

  if (length(curmchemdata$mz) < 1) stop("No mz data found!")

  result <- compute_chemical_score(
    mchemicaldata = curmchemdata,
    adduct_weights = adduct_weights,
    global_cor = global_cor,
    corthresh = corthresh,
    filter.by = filter.by,
    max_diff_rt = max_diff_rt,
    chemicalid = query$chemical_ID,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check
  )

  setwd("..")

  if (result$chemical_score >= (-100)) {
    result$filtdata <- result$filtdata[order(result$filtdata$mz), ]
    cur_chem_score <- result$chemical_score
    chemscoremat <- cbind(cur_chem_score, result$filtdata)
    chemscoremat <- as.data.frame(na.omit(chemscoremat))
  }
  return(chemscoremat)
}