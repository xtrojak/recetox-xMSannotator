test_that("Integration test: RT based clustering works", {
    peaks <- readRDS("test-data/qc_solvent.rda")
    peaks <- dplyr::rename(peaks, peak = feature)
    peaks <- dplyr::rename_with(
        peaks,
        ~ paste0("intensity_", .x),
        starts_with("Tribrid")
    )
    peaks$peak <- as.integer(peaks$peak)

    intensity_matrix <- t(select(peaks, starts_with("intensity")))
    intensity_matrix <- magrittr::set_colnames(intensity_matrix, peaks$peak)
    correlation_matrix <- WGCNA::cor(intensity_matrix, use = "p", method = "p")

    correlation_threshold <- 0.7
    deep_split <- 2
    min_cluster_size <- 10
    network_type <- "unsigned"

    peak_modules <- compute_peak_modules(
        peak_intensity_matrix = intensity_matrix,
        peak_correlation_matrix = correlation_matrix,
        correlation_threshold = correlation_threshold,
        deep_split = deep_split,
        min_cluster_size = min_cluster_size,
        network_type = network_type
    )

    peaks <- inner_join(peaks, peak_modules, by = "peak")
    rt_modules <- compute_rt_modules(peaks)
    expect_equal(nrow(peaks), nrow(rt_modules))
})