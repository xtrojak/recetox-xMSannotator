library(xmsannotator)
library(dplyr)
library(tibble)

compute_rt_modules <- function(peak_table, peak_modules) {
    peaks <- inner_join(peak_table, peak_modules, by = "peak")
    modules <- unique(peak_modules$module)
    rt_cluster <- {}
    for (id in modules) {
        subdata <- peaks %>% slice(filter(peak_modules, peak_modules$module == id)$peak)
        subdata$rt_cluster <- meanShiftR::meanShift(subdata$rt, iterations = 100)$assignment
        rt_cluster <- bind_rows(rt_cluster, subdata %>% select(peak, rt_cluster))
    }
    return(rt_cluster)
}

load("/home/hechth/dev/git/hechth/recetox-xMSannotator/xmsannotator/tests/testthat/test-data/sample_peaks.rda")
peaks$peak <- as.integer(rownames(peaks))
names(peaks) <- c("mz", "rt", "intensity_rep1", "intensity_rep2", "intensity_rep3", "peak")
peak_table <- peaks
peaks <- NA

peak_intensity_matrix <- t(select(peak_table, starts_with("intensity")))
peak_intensity_matrix <- magrittr::set_colnames(peak_intensity_matrix, peak_table$peak)
peak_correlation_matrix <- WGCNA::cor(peak_intensity_matrix, use = "p", method = "p")

correlation_threshold = 0.7
deep_split = 2
min_cluster_size = 10
network_type = "unsigned"

peak_modules <- compute_peak_modules(
  peak_intensity_matrix = peak_intensity_matrix,
  peak_correlation_matrix = peak_correlation_matrix,
  correlation_threshold = correlation_threshold,
  deep_split = deep_split,
  min_cluster_size = min_cluster_size,
  network_type = network_type
)

stage1 <- compute_rt_modules(peak_table, peak_modules)
names(stage1)