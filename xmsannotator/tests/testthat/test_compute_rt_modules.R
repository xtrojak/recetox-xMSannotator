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

test_that("scores are comaparble", {
    peaks <- readRDS("test-data/xmsannotator_qc_matrix_stage1.rda")
    peaks <- dplyr::rename_with(
        peaks,
        ~ paste0("intensity_", .x),
        starts_with("Tribrid")
    )
    peaks$peak <- as.integer(rownames(peaks))
    peaks <- dplyr::rename(
        peaks,
        rt = time,
        module = Module,
        RTclust_old = RTclust
    )

    matrix_with_clust <- compute_rt_modules(peaks, peak_width = 1)
    y2d <- gplots::hist2d(
        matrix_with_clust$RTClust_old,
        matrix_with_clust$RTclust
    )
    mutual_information <- entropy::mi.empirical(y2d$counts)
    expect_true(mutual_information > 2.5)
})
