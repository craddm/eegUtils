#' Plot ERP images
#'
#' Plot an ERP image from a single electrode. Smooths over a series of trials in
#' order to make across-trial patterns more apparent.
#'
#' @examples
#' erp_image(demo_epochs, electrode = "A31")
#' erp_image(demo_epochs, electrode = "A31", interpolate = TRUE)
#' erp_image(demo_epochs, electrode = "A31", smoothing = 5)
#' @param data Data frame to be plotted. Requires an amplitude column.
#' @param ... Other arguments passed to the method.
#' @import ggplot2
#' @importFrom scales squish
#' @export

erp_image <- function(data,
                      ...){
  UseMethod("erp_image", data)
}

#' @param electrode electrode at which to generate an ERP image.
#' @param smoothing Number of trials to smooth over when generating image
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#' @param interpolate Perform interpolation to produce smoother looking plots. Defaults to FALSE.
#' @describeIn erp_image Default function operates on normal data frames
#' @export
erp_image.default <- function(data,
                              electrode = "Cz",
                              smoothing = 10,
                              clim = NULL,
                              interpolate = FALSE,
                              ...) {

  required_cols <- c("electrode", "time", "amplitude", "epoch")
  col_names <- names(data)

  if (!all(required_cols %in% col_names)) {
    stop("Required columns ",
         required_cols[!required_cols %in% col_names], "missing.")
  }

  if (all(electrode %in% data$electrode)) {
    create_erpimage(data,
                    electrode = electrode,
                    smoothing = smoothing,
                    clim = clim)
  } else {
    stop("Electrode not found.")
  }

}

#'@describeIn erp_image Create an \code{erp_image} from \code{eeg_epochs}
#'@export
erp_image.eeg_epochs <- function(data,
                      electrode = "Cz",
                      smoothing = 10,
                      clim = NULL,
                      interpolate = FALSE,
                      ...) {

   if (!electrode %in% names(data$signals)) {
    stop("Specified electrode not found.")
  }
  data <- select_elecs(data,
                       electrode = electrode)
  data <- as.data.frame(data,
                          long = TRUE)
  create_erpimage(data,
                  electrode = electrode,
                  smoothing = smoothing,
                  clim = clim,
                  interpolate = interpolate)
}

#' @param component \code{eeg_ICA} component to plot
#' @describeIn erp_image Plot component image from \code{eeg_ICA}
#' @export
erp_image.eeg_ICA <- function(data,
                              component = "Comp1",
                              smoothing = 10,
                              clim = NULL,
                              interpolate = FALSE,
                              ...) {

  if (!component %in% names(data$signals)) {
    stop("Specified component not found.")
  }
  data <- select_elecs(data,
                       component = component)
  data <- as.data.frame(data,
                        long = TRUE)
  create_erpimage(data,
                  electrode = component,
                  smoothing = smoothing,
                  clim = clim,
                  interpolate = interpolate)
}

#' Function for creating an ERP image
#'
#' @param data Data frame to be plotted. Requires an amplitude column.
#' @param electrode electrode at which to generate an ERP image.
#' @param smoothing Number of trials to smooth over when generating image
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#' @param interpolate Turn on geom_raster() interpolation for smoother images.
#' @keywords internal
create_erpimage <- function(data,
                            electrode,
                            smoothing,
                            clim,
                            interpolate = FALSE) {

  n_times <- length(unique(data$time))
  n_epochs <- length(unique(data$epoch))
  sel_rows <- data$electrode %in% electrode
  data <- data[sel_rows, ]
  data$smooth_time <- rep(seq(min(data$time),
                              max(data$time),
                              length.out = n_times),
                          times = n_epochs)
  data$smooth_amp <- as.numeric(stats::filter(data$amplitude,
                                              rep(1 / smoothing,
                                                  smoothing),
                                              sides = 2))
  data$epoch <- as.numeric(factor(data$epoch))
  if (is.null(clim)) {
    clim <- max(abs(max(data$smooth_amp, na.rm = TRUE)),
                abs(min(data$smooth_amp, na.rm = TRUE)))
    clim <- c(-clim, clim)
  } else if (length(clim) != 2) {
    clim <- max(abs(max(data$smooth_amp, na.rm = TRUE)),
                abs(min(data$smooth_amp, na.rm = T)))
    clim <- c(-clim, clim)
  }

  ggplot2::ggplot(data,
                  aes(x = smooth_time,
                      y = epoch,
                      fill = smooth_amp)) +
    geom_raster(interpolate = interpolate) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               size = 1) +
    scale_fill_distiller(palette = "RdBu",
                         limits = clim,
                         oob = scales::squish) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    labs(x = "Time (s)", fill = "Amplitude", y = "Epoch number") +
    ggtitle(paste("ERP Image for electrode", electrode))
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
#' @examples
#' erp_raster(demo_epochs)
#' erp_raster(demo_epochs, interpolate = TRUE)
#' @param data An \code{eeg_epochs} object
#' @param anat_order Arrange the channels in a more anatomically representative
#'   order. Defaults to TRUE.
#' @param time_lim Time limits of plot - should be a character vector (e.g.
#'   c(-.2, .5))
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#' @param interpolate Smooth the raster plot. Defaults to FALSE.
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom scales squish
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export

erp_raster <- function(data,
                       anat_order = TRUE,
                       time_lim = NULL,
                       clim = NULL,
                       interpolate = FALSE) {

  if (!is.null(time_lim)){
    data <- select_times(data, time_lim)
  }

  if (all(anat_order && !is.null(data$chan_info))) {
    chan_order <- arrange_chans(data)
    data$signals <- data$signals[, chan_order]
  }

  data <- data.frame(data$signals,
                     time = data$timings$time)
  data <- split(data,
                data$time)
  data <- lapply(data,
                 Matrix::colMeans)
  data <- as.data.frame(do.call(rbind, data))
  data <- tidyr::gather(data,
                        electrode,
                        amplitude,
                        -time,
                        factor_key = TRUE)
  if (is.null(clim)) {
    clim <- c(min(data$amplitude),
              max(data$amplitude))
  }
  ggplot2::ggplot(data, aes(x = time,
                            y = electrode,
                            fill = amplitude)) +
    geom_raster(interpolate = interpolate) +
    geom_vline( xintercept = 0,
                linetype = "dashed",
                size = 2) +
    scale_fill_distiller(palette = "RdBu",
                         limits = c(clim[1],
                                    clim[2]),
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
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

arrange_chans <- function(data) {

  # Pick out electrodes on the midline (theta = 180 or 0)
  midline <- data$chan_info$pol_theta == 180 | data$chan_info$pol_theta == 0
  midline_labels <- data$chan_info$electrode[midline]
  midline_dist <- ifelse(data$chan_info$angle[midline] == 180,
                         -data$chan_info$radius[midline],
                         data$chan_info$radius[midline])

  # Order them from back to front
  midline_labels <- midline_labels[sort(midline_dist,
                                        decreasing = T,
                                        index.return = T)$ix]

  # Pick out electrodes from the left hemisphere
  left <- sapply(data$chan_info$pol_theta,
                 function(x) x < 0 && x > -180)
  left_labels <- data$chan_info$electrode[left]
  left_dist <- data$chan_info$y[left]
  left_labels <- left_labels[sort(left_dist,
                                  index.return = T)$ix]

  # Pirck out electrodes from the right hemisphere
  right <- sapply(data$chan_info$pol_theta,
                  function(x) x > 0 && x < 180)
  right_labels <- data$chan_info$electrode[right]
  right_dist <- data$chan_info$y[right]
  right_labels <- right_labels[sort(right_dist,
                                    index.return = T)$ix]

  all_labels <- c(left_labels,
                  midline_labels,
                  right_labels)
  elecs_found <- toupper(names(data$signals)) %in% all_labels
  elecs_not_found <- !elecs_found
  new_ord <- match(toupper(c(all_labels,
                             names(data$signals[!elecs_found]))),
                   toupper(names(data$signals)))
  new_ord
}
