#' Plot ERP images
#'
#' @param df data frame to be plotted.
#'
#' @import ggplot2
#' @import dplyr
#' @import viridis
#' @export
#'
#'

erp_image <- function(df, electrode = "Cz") {
  n_times <- length(unique(df$time))

  df <- df %>%
    filter(electrode == get("electrode")) %>%
    mutate(smooth_time = rep((seq(time[1], time[n_times], length.out = n_times)), times = length(unique(epoch))))

  ggplot(df,aes(x = smooth_time, y = epoch, fill = amplitude)) +
    geom_raster(interpolate = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
    scale_fill_viridis(limits = c(-20, 20)) +
    theme_bw() +
    xlab("Time (ms)")
}
