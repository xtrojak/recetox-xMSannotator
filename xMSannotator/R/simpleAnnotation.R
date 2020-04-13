simpleAnnotation <- function(dataA, max.mz.diff = 10, num_nodes = 2, queryadductlist = c("M+2H", "M+H+NH4", "M+ACN+2H", "M+2ACN+2H", 
    "M+H", "M+NH4", "M+Na", "M+ACN+H", "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H"), gradienttype = "Acetonitrile", mode = "pos", 
    outloc, db_name = "KEGG") {
    if (db_name == "KEGG") {
        data(keggCompMZ)
        chemCompMZ <- keggCompMZ
        suppressWarnings(rm(keggCompMZ))
    } else if (db_name == "HMDB") {
        data(hmdbCompMZ)
        chemCompMZ <- hmdbCompMZ
        suppressWarnings(rm(hmdbCompMZ))
    } else if (db_name == "T3DB") {
        data(t3dbCompMZ)
        chemCompMZ <- t3dbCompMZ
        suppressWarnings(rm(t3dbCompMZ))
    } else if (db_name == "LipidMaps") {
        data(lipidmapsCompMZ)
        chemCompMZ <- lipidmapsCompMZ
        suppressWarnings(rm(lipidmapsCompMZ))
    }

    data(adduct_table)
    adduct_table <- as.data.frame(adduct_table)
    adduct_table <- unique(adduct_table)
    
    WGCNA::allowWGCNAThreads(nThreads = num_nodes)

    suppressWarnings(dir.create(outloc))
    
    suppressWarnings(if (queryadductlist == "all" & mode == "pos") {
        adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "positive") | (adduct_table$Type == gradienttype & adduct_table$Mode == "positive")]
        adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
    } else if (queryadductlist == "all" & mode == "neg") {
        adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "negative") | (adduct_table$Type == gradienttype & adduct_table$Mode == "negative")]
        adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
    } else {
        adduct_names <- adduct_table$Adduct[which(adduct_table$Adduct %in% queryadductlist)]
        adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
        
        if (length(adduct_names) < 1) {
            stop("Invalid adducts selected.")
        }
    })

    adduct_names <- unique(adduct_names)
    chemCompMZ <- chemCompMZ[which(chemCompMZ$Adduct %in% adduct_names), ]
    
    print("Dimension of database matrix:")
    print(dim(chemCompMZ))

    cl <- parallel::makeCluster(num_nodes)

    parallel::clusterEvalQ(cl, library(XML))
    parallel::clusterEvalQ(cl, library(R2HTML))
    parallel::clusterEvalQ(cl, library(RCurl))
    parallel::clusterEvalQ(cl, library(SSOAP))
    parallel::clusterEvalQ(cl, library(limma))
    parallel::clusterEvalQ(cl, library(plyr))
    parallel::clusterEvalQ(cl, library(png))
    
    parallel::clusterEvalQ(cl, "SSOAP::processWSDL")

    parallel::clusterExport(cl, "find.Overlapping.mzs")
    parallel::clusterExport(cl, "find.Overlapping.mzsvparallel")
    parallel::clusterExport(cl, "overlapmzchild")
    parallel::clusterExport(cl, "getVenn")
    parallel::clusterExport(cl, "adduct_table")
    
    print("Mapping m/z to metabolites:")
    
    dataA <- dataA[, 1:2]
    adduct_names <- as.character(adduct_names)

    if (length(adduct_names) > 1) {
        l2 <- parallel::parLapply(cl, seq_along(adduct_names), xMSannotator::Annotationbychemical_IDschild, dataA = dataA, queryadductlist = adduct_names, adduct_type = c("S",
            gradienttype), max.mz.diff = max.mz.diff, outloc = outloc, keggCompMZ = chemCompMZ, otherdbs = FALSE, otherinfo = FALSE, 
            adduct_table = adduct_table, num_nodes = num_nodes)

        levelB_res <- {
        }
        for (j in seq_along(l2)) {
            if (length(l2[[j]]) > 1) {
                levelB_res <- rbind(levelB_res, l2[[j]])
            }
        }
        
        rm(l2)
    } else {
        levelB_res <- xMSannotator::Annotationbychemical_IDschild(adduct_index = 1, dataA = dataA, queryadductlist = adduct_names, adduct_type = c("S",
            gradienttype), max.mz.diff = max.mz.diff, outloc = outloc, keggCompMZ = chemCompMZ, otherdbs = FALSE, otherinfo = FALSE, 
            adduct_table = adduct_table, num_nodes = num_nodes)
    }

    parallel::stopCluster(cl)

    if (length(levelB_res) < 2) {
        print("No matches found.")
        return(-1)
    }
    
    levelB_res$mz <- as.numeric(as.character(levelB_res$mz))
    levelB_res$time <- as.numeric(as.character(levelB_res$time))
    levelB_res <- as.data.frame(levelB_res)

    uniq_formula <- as.character(unique(levelB_res$Formula))
    bad_formula <- which(is.na(uniq_formula) == TRUE)
    if (length(bad_formula) > 0) {
        uniq_formula <- uniq_formula[-bad_formula]
    }
    
    cl <- parallel::makeCluster(num_nodes)
    on.exit(parallel::stopCluster(cl))

    levelB_res_check <- parallel::parLapply(cl, uniq_formula, function(formula) {
        xMSannotator::check_golden_rules(curformula = as.character(formula), NOPS_check = TRUE)
    })
    
    levelB_res_check2 <- plyr::ldply(levelB_res_check, rbind)
    levelB_res_check3 <- levelB_res_check2[which(levelB_res_check2[, 2] == 1), ]
    levelB_res <- levelB_res[which(levelB_res$Formula %in% as.character(levelB_res_check3[, 1])), ]
    
    water_adducts <- c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
    water_adduct_ind <- which(levelB_res$Adduct %in% water_adducts)
    
    if (length(water_adduct_ind) > 0) {
        levelB_res2 <- levelB_res[water_adduct_ind, ]
        levelB_res <- levelB_res[-water_adduct_ind, ]
        
        sind1 <- seq(1:dim(levelB_res2)[1])
        
        levelB_res_check3 <- parallel::parLapply(cl, sind1, function(j) {
            curformula <- as.character(levelB_res2$Formula[j])
            numoxygens <- xMSannotator::check_element(curformula, "O")

            res <- cbind(curformula, numoxygens > 0)
            res <- as.data.frame(res)
        })
        
        levelB_res_check4 <- plyr::ldply(levelB_res_check3, rbind)

        if (length(which(levelB_res_check4[, 2] == 1)) > 0) {
            levelB_res_check4 <- levelB_res_check4[which(levelB_res_check4[, 2] == 1), ]
            valid_form <- which(levelB_res2$Formula %in% as.character(levelB_res_check4[, 1]))

            if (length(valid_form) > 0) {
                levelB_res <- rbind(levelB_res, levelB_res2[valid_form, ])
            }
        }
    }

    dupmz <- levelB_res$mz[which(duplicated(levelB_res$mz) == TRUE)]
    MatchCategory <- rep("Multiple", dim(levelB_res)[1])
    MatchCategory[-which(levelB_res$mz %in% dupmz)] <- "Unique"
    levelB_res <- cbind(MatchCategory, levelB_res)

    rownames(levelB_res) <- NULL
    return(levelB_res)
}
