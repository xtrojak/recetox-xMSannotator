#' @import dplyr
#' @export
compute_rt_modules <- function(peak_table, peak_modules, bandwidth = 1) {
    peaks <- inner_join(peak_table, peak_modules, by = "peak")
    modules <- unique(peak_modules$module)
    rt_cluster <- {}
    for (id in modules) {
        subdata <- peaks %>% slice(filter(peak_modules, peak_modules$module == id)$peak)
        subdata$rt_cluster <- meanShiftR::meanShift(
            subdata$rt,
            iterations = 100,
            bandwidth = bandwidth
        )$assignment
        rt_cluster <- bind_rows(rt_cluster, subdata %>% select(peak, rt_cluster))
    }
    return(rt_cluster)
}
