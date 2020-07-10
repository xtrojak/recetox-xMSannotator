getVenn <- function(dataA, name_a, dataB, name_b, mz.thresh = 10, time.thresh = 30, alignment.tool = NA, xMSanalyzer.outloc, use.unique.mz = FALSE, plotvenn = TRUE, num_nodes = 2) {
  dir.create(xMSanalyzer.outloc, showWarnings = FALSE)

  plotvenn <- FALSE
  data_a <- as.data.frame(dataA)
  data_b <- as.data.frame(dataB)
  rm(dataA)
  rm(dataB)

  if (use.unique.mz == TRUE) {
    data_a <- find.Unique.mzs.sameset(dataA = data_a, dataB = data_a, mz.thresh = mz.thresh, time.thresh = time.thresh, alignment.tool = alignment.tool)
    data_a <- data_a$uniqueA

    data_b <- find.Unique.mzs.sameset(dataA = data_b, dataB = data_b, mz.thresh = mz.thresh, time.thresh = time.thresh, alignment.tool = alignment.tool)
    data_b <- data_b$uniqueA
  }

  common <- find.Overlapping.mzsvparallel(data_a, data_b, mz.thresh, time.thresh = time.thresh, alignment.tool = alignment.tool, num_nodes = num_nodes)

  if (is.na(time.thresh) == FALSE) {
    mznames <- c("index.A", "mz.data.A", "time.data.A", "index.B", "mz.data.B", "time.data.B", "time.difference")
  }else {
    mznames <- c("index.A", "mz.data.A", "index.B", "mz.data.B")
  }

  if (length(common) > 0) {
    colnames(common) <- mznames
  }

  return(list("common" = common))
}
