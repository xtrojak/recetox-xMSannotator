#' Plot the clustering and density function.
#'
#' @param pdf Density estimate using stats::density(...) function. Black line.
#' @param data Raw data used to estimate kernel density.
#' @param peak_indices Computed indices where to plot the peaks.
#' Red circles & lines.
#'
#' @return
#' @export
#'
#' @examples
plot_clustering <- function(pdf, data, peak_indices, plot_lines = FALSE) {
    plot(pdf, type = "l", xaxt = "n")
    axis(1, at = data, las = 2)
    points(pdf$x[peak_indices], pdf$y[peak_indices], col = "red")

    if (plot_lines) abline(v = pdf$x[peak_indices], col = "red")
}


#' Compute the positions of peak indices for given kernel density estimate.
#'
#' @param density Kernel density estimates (stats::density$y).
#'
#' @return Indices at which kernel density estimate position is a cluster.
#' @export
#' @import pastecs
#' @import rlist
#'
#' @examples
compute_peak_indices <- function(density) {
    turning_points <- turnpoints(ts(density), calc.proba = FALSE)
    peak_indices <- turning_points$pos[turning_points$peaks]
    filtered <- list.filter(peak_indices, density[.] > .Machine$double.eps)
    return(filtered)
}


#' Compute number of equally spaced kernel points for density estimation.
#'
#' @param data Data for which to estimate the number of kernel points.
#'
#' @return Nearest power of 2 to length(data), but at least 2048.
#' @export
#'
#' @examples
compute_num_kernel_points <- function(data) {
    num_points <- length(data)
    nearest_power_of_two <- 2**ceiling(log2(num_points))
    max(nearest_power_of_two, 2048)
}


#' Compute kernel density estimate.
#'
#' @param data Data for which to estimate the kernel density.
#' @param width Kernel width to use. Default 1. [0.5;2]
#' @param kernel Kernel to use. See stats::density for details.
#' Recommended ["gaussian","biweight","cosine","epanechnikov"].
#'
#' @return Kernel density estimate (stats::density).
#' @export
#'
#' @examples
estimate_kernel_density <- function(data, width = 1, kernel = "gaussian") {
    start <- max(min(data) - 1, 0)
    end <- max(data) + 1

    pdf <- density(
        data,
        from = start,
        to = end,
        n = compute_num_kernel_points(data),
        width = width,
        kernel = kernel
    )
}

#' Compute cluster positions in data
#'
#' @param data Data for which to compute the clustering.
#' @param width Bandwidth of the clustering.
#' Scan rate should be consulted to optimize clustering. [0.5;2.0]
#'
#' @return Positions of dense clusters.
#' @export
#'
#' @examples
compute_cluster_positions <- function(data, width = 1, kernel = "gaussian", do_plot = FALSE) {
    pdf <- estimate_kernel_density(data, width = width, kernel = kernel)
    peak_indices <- compute_peak_indices(pdf$y)
    if (do_plot) plot_clustering(pdf, data, peak_indices)

    peak_positions <- pdf$x[peak_indices]
    return(peak_positions)
}

#' Compute the assignment of cluster ID based on nn-algorithm.
#'
#' @param clusters Position of cluster centers.
#' @param data Data to assign to clusters
#'
#' @return Indices of closest cluster for each query point
#' @import RANN
#' @export
#'
#' @examples
compute_cluster_assignments <- function(clusters, data) {
    nn2(clusters, data, k = 1)$nn.idx[, 1]
}

#' Compute RT modules based on intensity modules and chromatographic peak width.
#'
#' @param peak_table Feature table with columns ['rt', 'module','peak'].
#' @param peak_width Estimated chromatographic peak width.
#'
#' @return
#' @export
#' @import dplyr
#'
#' @examples
compute_rt_modules <- function(peak_table, peak_width = 1) {
    modules <- unique(peak_table$module)
    rt_cluster <- {}

    for (id in modules) {
        subdata <- peak_table %>% filter(module == id)
        cluster_positions <- compute_cluster_positions(subdata$rt, width = peak_width, do_plot = TRUE)
        subdata$RTclust <- compute_cluster_assignments(cluster_positions, subdata$rt)

        rt_cluster <- bind_rows(rt_cluster, subdata %>% select(peak, RTclust))
    }

    return(left_join(peak_table, rt_cluster, by = "peak"))
}