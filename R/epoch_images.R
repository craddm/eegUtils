#' Plot ERP images
#'
#' @param data Data frame to be plotted. Requires an amplitude column.
#' @param electrode electrode at which to generate an ERP image.
#' @param smoothing Number of trials to smooth over when generating image
#' @param clim Character vector of min and max values of plotting colour range. e.g. c(-5,5). Defaults to min and max.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @import ggplot2
#' @import dplyr
#' @import viridis
#' @export
#'
#'

erp_image <- function(data, electrode = "Cz", smoothing = 10, clim = NULL) {

  n_times <- length(unique(data$time))

  data <- filter(data, electrode == get("electrode"))
  data <- mutate(data,
                 smooth_time = rep((seq(time[1], time[n_times], length.out = n_times)), times = length(unique(epoch))),
                 smooth_amp = as.numeric(stats::filter(amplitude, rep(1/smoothing,smoothing), sides = 2)),
                 epoch = as.numeric(factor(epoch)))

  if (is.null(clim)) {
    clim <- max(abs(max(data$smooth_amp, na.rm = T)), abs(min(data$smooth_amp, na.rm = T))) %>% c(-., .)
  } else if (length(clim) != 2) {
    clim <- max(abs(max(data$smooth_amp, na.rm = T)), abs(min(data$smooth_amp, na.rm = T))) %>% c(-., .)
  }

  ggplot(data,aes(x = smooth_time, y = epoch, fill = smooth_amp)) +
    geom_raster(interpolate = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
    scale_fill_distiller(palette = "RdBu", limits = c(clim[1],clim[2]), oob = scales::squish) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_classic() +
    xlab("Time (ms)")
}
