test_that("RT based clustering works", {

    load("test-data/sample_peaks.rda")
    peaks$peak <- as.integer(rownames(peaks))
    names(peaks) <- c("mz", "rt", "intensity_rep1", "intensity_rep2", "intensity_rep3", "peak")
    peak_table <- peaks
    peaks <- NA

    peak_intensity_matrix <- t(select(peak_table, starts_with("intensity")))
    peak_intensity_matrix <- magrittr::set_colnames(peak_intensity_matrix, peak_table$peak)
    peak_correlation_matrix <- WGCNA::cor(peak_intensity_matrix, use = "p", method = "p")

    correlation_threshold <- 0.7
    deep_split <- 2
    min_cluster_size <- 10
    network_type <- "unsigned"

    peak_modules <- compute_peak_modules(
        peak_intensity_matrix = peak_intensity_matrix,
        peak_correlation_matrix = peak_correlation_matrix,
        correlation_threshold = correlation_threshold,
        deep_split = deep_split,
        min_cluster_size = min_cluster_size,
        network_type = network_type
    )

    rt_modules <- compute_rt_modules(peak_table, peak_modules, bandwidth = 10)

    peaks <- left_join(peak_table, peak_modules, by = "peak")
    stage1 <- left_join(peaks, rt_modules, by = "peak")
    stage1$Module_RTclust <- stringr::str_c(stage1$module, stage1$rt_cluster, sep = "_")
    write.csv(stage1, file = "stage1_rcx_final.csv")
})