overlapmzchild <- function(j, mz.thresh = 10, time.thresh = NA, data_a, data_b) {
  commat <- { }
  ppmb <- (mz.thresh) * (data_a$mz[j] / 1000000)

  getbind_same <- which(abs(data_b$mz - data_a$mz[j]) <= ppmb)

  if (is.na(time.thresh) == FALSE) {
    if (length(getbind_same) > 0) {
      nearest_time_diff <- 10000
      bestmatch <- { }
      rnames <- { }
      temp <- { }
      commat <- { }

      commat <- lapply(seq_along(getbind_same), function(comindex) {
        tempA <- cbind(j, data_a[j, c(1, 2)])
        tempB <- cbind(getbind_same[comindex], data_b[getbind_same[comindex], c(1, 2)])
        temp <- cbind(tempA, tempB)

        timediff <- abs(data_a[j, 2] - data_b[getbind_same[comindex], 2])

        temp <- cbind(temp, timediff)

        if (timediff < time.thresh && timediff <= nearest_time_diff) {
          bestmatch <- as.data.frame(temp)
          nearest_time_diff <- timediff
          temp <- as.data.frame(temp)
        }
        return(temp)
      })

      commat <- plyr::ldply(commat, rbind)
    }
  } else if (length(getbind_same) > 0) {
    commat <- lapply(seq_along(getbind_same), function(comindex) {
      cbind(j, data_a[j, 1], getbind_same[comindex], data_b[getbind_same[comindex], 1])
    })
    commat <- plyr::ldply(commat, rbind)
    commat <- as.data.frame(commat)
  }

  return(as.data.frame(commat))
}
