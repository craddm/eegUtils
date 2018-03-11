#' Plot ERP images
#'
#' Plot an ERP image from a single electrode. Smooths over a series of trials in
#' order to make across-trial patterns more apparent.
#'
#' @param data Data frame to be plotted. Requires an amplitude column.
#' @param electrode electrode at which to generate an ERP image.
#' @param smoothing Number of trials to smooth over when generating image
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @import ggplot2
#' @import dplyr
#' @importFrom scales squish
#' @export
#'


erp_image <- function(data, electrode = "Cz", smoothing = 10, clim = NULL) {

  if (is.eeg_data(data)) {
    data <- as.data.frame(data, long = TRUE)
  }

  n_times <- length(unique(data$time))

  data <- dplyr::filter(data, electrode == !!electrode)
  data <- dplyr::mutate(data,
                 smooth_time = rep((seq(time[1],
                                        time[n_times],
                                        length.out = n_times)),
                                   times = length(unique(epoch))),
                 smooth_amp = as.numeric(stats::filter(amplitude,
                                                       rep(1 / smoothing,
                                                           smoothing),
                                                       sides = 2)),
                 epoch = as.numeric(factor(epoch)))

  if (is.null(clim)) {
    clim <- max(abs(max(data$smooth_amp, na.rm = T)),
                abs(min(data$smooth_amp, na.rm = T))) %>% c(-., .)
  } else if (length(clim) != 2) {
    clim <- max(abs(max(data$smooth_amp, na.rm = T)),
                abs(min(data$smooth_amp, na.rm = T))) %>% c(-., .)
  }

  ggplot2::ggplot(data, aes(x = smooth_time, y = epoch, fill = smooth_amp)) +
    geom_raster(interpolate = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
    scale_fill_distiller(palette = "RdBu",
                         limits = c(clim[1], clim[2]),
                         oob = scales::squish) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    xlab("Time (s)")
}

#' ERP raster plot
#'
#' Create a plot showing the ERP at every channel as a single ERP image. By
#' default, this attempts to find channel locations and rearrange the channels
#' such that spatial patterns on the scalp are more discernible. It orders the
#' rows from the most posterior electrode on the right hemisphere to the most
#' anterior electrode on the left hemisphere, with midline electrodes in the
#' middle. If no locations are found, it simply displays the data in its
#' original order.
#'
#' @param data An \code{eeg_epochs} object
#' @param anat_order Arrange the channels in a more anatomically representative
#'   order.
#' @param time_lim Time limits of plot - should be a character vector (e.g.
#'   c(-.2, .5))
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#' @import ggplot2
#' @import dplyr
#' @importFrom scales squish
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @export

erp_raster <- function(data, anat_order = TRUE, time_lim = NULL, clim = NULL) {

  if (!is.null(time_lim)){
    data <- select_times(data, time_lim)
  }

  if (all(anat_order && !is.null(data$chan_info))) {
    chan_order <- arrange_chans(data)
    data$signals <- data$signals[, chan_order]
  }

  data <- data.frame(data$signals, time = data$timings$time)
  data <- dplyr::group_by(data, time)
  data <- dplyr::summarise_all(data, mean)
  data <- tidyr::gather(data, electrode, amplitude, -time, factor_key = TRUE)

  if (is.null(clim)) {
    clim <- c(min(data$amplitude), max(data$amplitude))
  }

  ggplot2::ggplot(data, aes(x = time, y = electrode, fill = amplitude)) +
    geom_raster(interpolate = TRUE) +
    geom_vline( xintercept = 0, linetype = "dashed", size = 2) +
    scale_fill_distiller(palette = "RdBu",
                         limits = c(clim[1], clim[2]),
                         oob = scales::squish) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    labs(x = "Time (s)")
}


#' Topographical channel ordering
#'
#' Rearrange channels in a more anatomical order (e.g. left to right, anterior
#' to posterior).
#'
#' @param data An \code{eeg_data} object
#' @return Vector of channel positions
#' @author Matt Craddock, \email{matt@mattcraddock.com}

arrange_chans <- function(data) {

  # Pick out electrodes on the midline (theta = 180 or 0)
  midline <- data$chan_info$theta == 180 | data$chan_info$theta == 0
  midline_labels <- data$chan_info$electrode[midline]
  midline_dist <- ifelse(data$chan_info$theta[midline] == 180,
                         -data$chan_info$radius[midline],
                         data$chan_info$radius[midline])

  # Order them from back to front
  midline_labels <- midline_labels[sort(midline_dist,
                                        decreasing = T,
                                        index.return = T)$ix]

  # Pick out electrodes from the left hemisphere
  left <- sapply(data$chan_info$theta, function(x) x < 0 && x > -180)
  left_labels <- data$chan_info$electrode[left]
  left_dist <- data$chan_info$y[left]
  left_labels <- left_labels[sort(left_dist,
                                  index.return = T)$ix]

  # Pirck out electrodes from the right hemisphere
  right <- sapply(data$chan_info$theta, function(x) x > 0 && x < 180)
  right_labels <- data$chan_info$electrode[right]
  right_dist <- data$chan_info$y[right]
  right_labels <- right_labels[sort(right_dist,
                                    index.return = T)$ix]

  all_labels <- c(left_labels, midline_labels, right_labels)
  elecs_found <- toupper(names(data$signals)) %in% all_labels
  elecs_not_found <- !elecs_found
  new_ord <- match(toupper(c(all_labels,
                             names(data$signals[!elecs_found]))),
                   toupper(names(data$signals)))
  new_ord
}
