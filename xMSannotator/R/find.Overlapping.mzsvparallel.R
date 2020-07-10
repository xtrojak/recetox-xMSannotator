find.Overlapping.mzsvparallel <-
  function(dataA, dataB, mz.thresh = 10, time.thresh = NA, alignment.tool = NA, num_nodes = 2) {
    data_a <- as.data.frame(dataA)
    data_b <- as.data.frame(dataB)
    rm(dataA)
    rm(dataB)

    col.names.dataA = colnames(data_a)
    col.names.dataB = colnames(data_b)

    if (alignment.tool == "apLCMS") {
      sample.col.start = 5
    } else if (alignment.tool == "XCMS") {
      sample.col.start = 9
      col.names.dataA[1] = "mz"
      col.names.dataA[2] = "time"
      col.names.dataB[1] = "mz"
      col.names.dataB[2] = "time"
      colnames(data_a) = col.names.dataA
      colnames(data_b) = col.names.dataB
    } else if (is.na(alignment.tool)) {
      col.names.dataA[1] = "mz"
      col.names.dataB[1] = "mz"

      if (is.na(time.thresh) == FALSE) {
        col.names.dataA[2] = "time"
        col.names.dataB[2] = "time"
      }

      colnames(data_a) = col.names.dataA
      colnames(data_b) = col.names.dataB
    }

    data_a <- as.data.frame(data_a)
    data_b <- as.data.frame(data_b)
    colnames(data_a) = col.names.dataA
    colnames(data_b) = col.names.dataB

    #create header for the matrix with common features
    if (is.na(time.thresh) == FALSE) {
      mznames <- c("index.A", "mz.data.A", "time.data.A", "index.B", "mz.data.B", "time.data.B", "time.difference")
    } else {
      mznames <- c("index.A", "mz.data.A", "index.B", "mz.data.B")
    }

    #Step 1 Group features by m/zdim(data_a)[1]
    parallel_list_vec <- 1:dim(data_a)[1]

    cl <- parallel::makeSOCKcluster(num_nodes)
    parallel::clusterExport(cl, "ldply")

    mz_groups <- parallel::parLapply(cl, parallel_list_vec, overlapmzchild, mz.thresh = mz.thresh, time.thresh = time.thresh, data_a = data_a, data_b = data_b)

    parallel::stopCluster(cl)
    #Step 2 Sub-group features from Step 1 by Retention time
    #find the features with RT values within the defined range as compared to the query feature

    commat <- data.frame()

    if (length(mz_groups) > 0) {
      commat <- plyr::ldply(mz_groups, rbind)
    }else {
      stop("No mz_groups defined.")
    }

    if (is.null(dim(commat)) == FALSE) {
      commat <- as.data.frame(commat)
    }
    rm(data_a)
    rm(data_b)

    return(commat)
  }
