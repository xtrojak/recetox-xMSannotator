Annotationbychemical_IDs <- function(dataA, queryadductlist = c("M+H"), adduct_type = c("S", "Acetonitrile"), adduct_table, max.mz.diff = 10, outloc, numnodes = 10, otherdbs = FALSE, otherinfo = FALSE, keggCompMZ, syssleep = 0.01) {
    adduct_names <- as.character(adduct_table[, 1])
    adductlist <- adduct_table[, 4]
    mult_charge <- adduct_table[, 3]
    num_mol <- adduct_table[, 2]
    names(adductlist) <- as.character(adduct_names)
    names(mult_charge) <- as.character(adduct_names)
    names(num_mol) <- as.character(adduct_names)
    
    suppressWarnings(dir.create(outloc))
    
    setwd(outloc)
    
    cl <- parallel::makeCluster(numnodes)
    
    parallel::clusterEvalQ(cl, library(XML))
    parallel::clusterEvalQ(cl, library(R2HTML))
    parallel::clusterEvalQ(cl, library(RCurl))
    parallel::clusterEvalQ(cl, library(SSOAP))
    parallel::clusterEvalQ(cl, library(png))
    
    parallel::clusterEvalQ(cl, "processWSDL")
    parallel::clusterEvalQ(cl, "Annotationbychemical_IDschild")
    parallel::clusterEvalQ(cl, "keggCompMZ")
    
    parallel::clusterExport(cl, "find.Overlapping.mzs")
    parallel::clusterExport(cl, "find.Overlapping.mzsvparallel")
    parallel::clusterExport(cl, "overlapmzchild")
    parallel::clusterExport(cl, "getVenn")
    
    parallel::clusterEvalQ(cl, library(limma))
    
    mzlist <- dataA[, 1]
    mz_group <- 10
    s1 <- seq(1, length(queryadductlist))
    
    mz.annot.res <- parallel::parLapply(cl, s1, Annotationbychemical_IDschild, dataA, queryadductlist, adduct_type, adduct_table, max.mz.diff, outloc, otherdbs, otherinfo, keggCompMZ)
    if (is(mz.annot.res, "try-error")) {
        Sys.sleep(10)
        mz.annot.res <- parallel::parLapply(cl, s1, Annotationbychemical_IDschild, dataA, queryadductlist, adduct_type, adduct_table, max.mz.diff, outloc, otherdbs, otherinfo, keggCompMZ)
    }
    
    if (FALSE) {
        for (mind in seq(1, length(mzlist), mz_group)) {
            
            stopind <- mind + mz_group - 1
            if (stopind > length(mzlist)) {
                stopind <- length(mzlist)
            }
            s1 <- dataA[mind:stopind, ]
            s1 <- unique(s1)
            num_mz <- dim(s1)[1]
            if (num_mz%%50 > 0) {
                Sys.sleep(syssleep)
            } else {
                Sys.sleep(syssleep)
            }
            
            if (num_mz > 100) {
                repeat {
                  cur.annot.res <- parallel::parLapply(cl, s1, xMSannotator::Annotationbychemical_IDschild, queryadductlist, adduct_type, adduct_table, max.mz.diff, outloc, otherdbs, otherinfo, keggCompMZ)
                  if (is(cur.annot.res, "try-error")) {
                    Sys.sleep(10)
                    cur.annot.res <- parallel::parLapply(cl, s1, xMSannotator::Annotationbychemical_IDschild, queryadductlist, adduct_type, adduct_table, max.mz.diff, outloc, otherdbs, otherinfo, keggCompMZ)
                  } else {
                    break
                  }
                }
                
                mz.annot.res <- c(mz.annot.res, cur.annot.res)
            } else {
                rescur <- xMSannotator::Annotationbychemical_IDschild(s1, queryadductlist, adduct_type, adduct_table, max.mz.diff, outloc, otherdbs, otherinfo, keggCompMZ)
                if (length(rescur) > 0) {
                  rescur <- as.matrix(rescur)
                  if (dim(rescur)[2] == 1) {
                    rescur <- t(rescur)
                    rescur <- as.data.frame(rescur)
                  }
                  rescur <- as.data.frame(rescur)
                  mz.annot.res[[length(mz.annot.res) + 1L]] <- rescur
                }
            }
            if (mind%%10 > 0) {
                Sys.sleep(syssleep/2)
            } else {
                Sys.sleep(1)
                parallel::stopCluster(cl)
                cl <- parallel::makeCluster(numnodes)
                
                
                parallel::clusterEvalQ(cl, library(XML))
                parallel::clusterEvalQ(cl, library(RCurl))
                parallel::clusterEvalQ(cl, library(SSOAP))
                parallel::clusterEvalQ(cl, library(png))
                parallel::clusterEvalQ(cl, "processWSDL")
                parallel::clusterEvalQ(cl, "Annotationbychemical_IDschild")
            }
        }
    }
    parallel::stopCluster(cl)
    
    res <- {
    }
    
    if (length(mz.annot.res) > 0) {
        for (mzl in seq_along(mz.annot.res)) {
            res <- rbind(res, mz.annot.res[[mzl]])
        }
    }
    res <- unique(res)
    
    return(res)
}
