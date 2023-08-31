#' Plot ERP images
#'
#' Plot an ERP image from a single electrode. Uses a boxcar smooth over a series of trials in
#' order to make across-trial patterns more apparent.
#'
#' @examples
#' erp_image(demo_epochs, electrode = "A31")
#' erp_image(demo_epochs, electrode = "A31", interpolate = TRUE)
#' erp_image(demo_epochs, electrode = "A31", smoothing = 5)
#' erp_image(compute_tfr(demo_epochs,
#'  foi = c(4, 30), n_cycles = 3, n_freq = 20, verbose = FALSE, keep_trials = TRUE),
#'  electrode = "A31", freq_range = c(8, 12))
#' @param data Data frame or `eegUtils` object to be plotted.
#' @param ... Other arguments passed to the method.
#' @author Matt craddock \email{matt@@mattcraddock.com}
#' @import ggplot2
#' @return A `ggplot` object
#' @export

erp_image <- function(data,
                      ...) {
  UseMethod("erp_image", data)
}

#' @export
erp_image.default <- function(data,
                              ...) {
  stop("Not implemented for objects of class ",
       class(data))
}

#' @export
erp_image.eeg_evoked <- function(data,
                                 ...) {
  stop("`erp_image()` requires single-trial data, but you have passed an `eeg_evoked` object.")
}

#' @param electrode Electrode for which to generate an ERP image.
#' @param time_lim Time limits of plot.
#' @param smoothing Number of trials to smooth over when generating image.
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#' @param interpolate Perform interpolation to produce smoother looking plots.
#'   Defaults to FALSE.
#' @param na.rm Remove trials with NA amplitudes after smoothing. Defaults to TRUE.
#' @describeIn erp_image Default function operates on normal data frames
#' @export
erp_image.data.frame <- function(data,
                                 electrode = "Cz",
                                 time_lim = NULL,
                                 smoothing = 10,
                                 clim = NULL,
                                 interpolate = FALSE,
                                 na.rm = TRUE,
                                 ...) {

  required_cols <- c("electrode",
                     "time",
                     "amplitude",
                     "epoch")

  if (length(electrode) > 1) {
    stop("Currently, only one electrode can be plotted at a time.")
  }
  col_names <- names(data)

  if (!all(required_cols %in% col_names)) {
    stop("Required columns ",
         required_cols[!required_cols %in% col_names], "missing.")
  }

  if (!is.null(time_lim)) {
    data <-
      dplyr::filter(
        data,
        time > time_lim[1],
        time < time_lim[2]
        )
  }

  if (all(electrode %in% data$electrode)) {
    sel_elec <- electrode
    data <-
      dplyr::filter(
        data,
        electrode %in% sel_elec
      )
    create_erpimage(data,
                    electrode = electrode,
                    smoothing = smoothing,
                    clim = clim,
                    na.rm = na.rm)
  } else {
    stop("Electrode not found.")
  }

}



#'@describeIn erp_image Create an `erp_image` from `eeg_epochs`
#'@export
erp_image.eeg_epochs <- function(data,
                                 electrode = "Cz",
                                 time_lim = NULL,
                                 smoothing = 10,
                                 clim = NULL,
                                 interpolate = FALSE,
                                 na.rm = TRUE,
                                 ...) {
  if (length(electrode) > 1) {
    stop("Currently, only one electrode can be plotted at a time.")
  }

  if (!electrode %in% names(data$signals)) {
    stop("Specified electrode not found.")
  }

  if (!is.null(time_lim)) {
    data <- filter(data,
                   time > time_lim[1],
                   time < time_lim[2])
  }

  data <- select_elecs(data,
                       electrode = electrode)
  data <- as.data.frame(data,
                        long = TRUE,
                        coords = FALSE)

  create_erpimage(data,
                  electrode = electrode,
                  smoothing = smoothing,
                  clim = clim,
                  interpolate = interpolate,
                  na.rm = na.rm)
}

#' @param component `eeg_ICA` component to plot
#' @describeIn erp_image Plot component image from `eeg_ICA`
#' @export
erp_image.eeg_ICA <- function(data,
                              component = "Comp001",
                              smoothing = 10,
                              clim = NULL,
                              interpolate = FALSE,
                              na.rm = TRUE,
                              ...) {

  if (length(component) > 1) {
    stop("Only one `component` should be supplied")
  }
  data <- select_elecs(data,
                       component = component)
  data <- as.data.frame(data,
                        long = TRUE,
                        coords = FALSE)
  create_erpimage(data,
                  electrode = component,
                  smoothing = smoothing,
                  clim = clim,
                  interpolate = interpolate,
                  na.rm = na.rm)

}

#' @param freq_range A numeric vector specify the range of frequencies to
#'   average over. A single number will find the closest matching frequency. A
#'   vector of length two will match average over frequencies within that range.
#' @describeIn erp_image Plot component image from `eeg_tfr`
#' @export
erp_image.eeg_tfr <- function(data,
                              electrode = "Cz",
                              time_lim = NULL,
                              smoothing = 10,
                              clim = NULL,
                              interpolate = FALSE,
                              freq_range = NULL,
                              na.rm = TRUE,
                              ...) {

  if (length(unique(epochs(data)$epoch)) == 1) {
    stop("`erp_image()` requires an `eeg_tfr` object with more than one trial.")
  }

  if (length(electrode) > 1) {
    stop("Currently, only one electrode can be plotted at a time.")
  }

  data <- select_elecs(data,
                       electrode)

  if (!is.null(freq_range)) {
    data <- select_freqs(data,
                         freq_range)
  }

  if (!is.null(time_lim)) {
    data <- filter(data,
                   time > time_lim[1],
                   time < time_lim[2])
  }

  data <- as.data.frame(data,
                        long = TRUE,
                        coords = FALSE)



  create_tfrimage(data,
                  electrode = electrode,
                  smoothing = smoothing,
                  clim = clim,
                  interpolate = interpolate,
                  na.rm = na.rm)

}

#' Function for creating an ERP image
#'
#' @param data Data frame to be plotted. Requires an amplitude column.
#' @param electrode electrode at which to generate an ERP image.
#' @param smoothing Number of trials to smooth over when generating image
#' @param clim Character vector of min and max values of plotting colour range.
#'   e.g. c(-5,5). Defaults to min and max.
#' @param interpolate Turn on `geom_raster()` interpolation for smoother images.
#' @param na.rm Remove trials with NA amplitudes after smoothing. Defaults to TRUE.
#' @keywords internal
create_erpimage <- function(data,
                            electrode,
                            smoothing,
                            clim,
                            interpolate = FALSE,
                            na.rm) {

  n_times <- length(unique(data$time))
  n_epochs <- length(unique(data$epoch))

  data <- split(data,
                data$time)

  data <- lapply(data,
                      function(x) {
                        x$smooth_amp <-
                          as.numeric(stats::filter(x$amplitude,
                                                   rep(1 / smoothing,
                                                       smoothing),
                                                   sides = 2))
                        x
                      })

  data <- data.table::rbindlist(data)
  data$epoch <- as.numeric(factor(data$epoch))

  if (is.null(clim)) {
    clim <- max(abs(max(data$smooth_amp, na.rm = TRUE)),
                abs(min(data$smooth_amp, na.rm = TRUE)))
    clim <- c(-clim, clim)
  } else if (length(clim) != 2) {
    clim <- max(abs(max(data$smooth_amp, na.rm = TRUE)),
                abs(min(data$smooth_amp, na.rm = TRUE)))
    clim <- c(-clim, clim)
  }

  if (na.rm) {
    data <- data[!is.na(data$smooth_amp), ]
  }

  ggplot2::ggplot(tibble::as_tibble(data),
                  aes(x = time,
                      y = epoch,
                      fill = smooth_amp)) +
    geom_raster(interpolate = interpolate) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 1) +
    scale_fill_distiller(palette = "RdBu",
                         limits = clim,
                         oob = scales::squish) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    labs(x = "Time (s)", fill = "Amplitude", y = "Epoch number") +
    ggtitle(paste("ERP Image for", electrode))
}


create_tfrimage <- function(data,
                            electrode,
                            smoothing,
                            clim,
                            interpolate = FALSE,
                            freq_range = NULL,
                            na.rm) {

  data <- dplyr::group_by(data,
                          time, electrode, epoch)
  data <- dplyr::summarise(data,
                           power = mean(power, na.rm = TRUE))
  n_times <- length(unique(data$time))
  n_epochs <- length(unique(data$epoch))

  data <- split(data,
                data$time)

  data <- lapply(data,
                 function(x) {
                   x$smooth_amp <-
                     as.numeric(stats::filter(x$power,
                                              rep(1 / smoothing,
                                                  smoothing),
                                              sides = 2))
                   x
                 })

  data <- data.table::rbindlist(data)
  data$epoch <- as.numeric(factor(data$epoch))

  if (na.rm) {
    data <- data[!is.na(data$smooth_amp), ]
  }

  if (is.null(clim)) {
    clim <- max(abs(max(data$smooth_amp, na.rm = TRUE)),
                abs(min(data$smooth_amp, na.rm = TRUE)))
    clim <- c(0, clim)
  } else if (length(clim) != 2) {
    clim <- max(abs(max(data$smooth_amp, na.rm = TRUE)),
                abs(min(data$smooth_amp, na.rm = TRUE)))
    clim <- c(0, clim)
  }

  ggplot2::ggplot(tibble::as_tibble(data),
                  aes(x = time,
                      y = epoch,
                      fill = smooth_amp)) +
    geom_raster(interpolate = interpolate) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 1) +
    scale_fill_viridis_c(limits = clim,
                         oob = scales::squish) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    labs(x = "Time (s)",
         fill = "Power (a.u.)",
         y = "Epoch number") +
    ggtitle(paste("TFR Image for electrode",
                  electrode))
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
#' library(ggplot2)
#' erp_raster(demo_epochs)
#' erp_raster(demo_epochs, interpolate = TRUE)
#' erp_raster(rm_baseline(demo_epochs, c(-.1, 0)), interpolate = TRUE)
#' erp_raster(demo_spatial) + facet_wrap(~epoch_labels)
#' @param data An `eeg_epochs` object
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
#' @return A `ggplot` object
#' @export

erp_raster <- function(data,
                       anat_order = TRUE,
                       time_lim = NULL,
                       clim = NULL,
                       interpolate = FALSE) {

  if (
    inherits(
      data, c("eeg_tfr", "eeg_ICA")
      )
    ) {
    stop(
      paste(
        "Not currently implemented for objects of class",
        class(data)[[1]]
        )
      )
  }

  if (!is.null(time_lim)){
    data <- select_times(data, time_lim)
  }

  if (all(anat_order && !is.null(data$chan_info))) {
    chan_order <- arrange_chans(channels(data),
                                channel_names(data))
    chan_lab_order <- channel_names(data)[chan_order]
  } else {
    chan_lab_order <- channel_names(data)
  }

  data <- as.data.frame(data,
                        long = TRUE)
  data$electrode <- factor(data$electrode,
                           levels = chan_lab_order)
  if (is.null(clim)) {
    clim_max <- max(abs(range(data$amplitude)))
    clim <- c(-clim_max, clim_max) / 10
  }

  ggplot2::ggplot(data,
                  aes(x = time,
                      y = electrode,
                      fill = amplitude)) +
    stat_summary_by_fill(interpolate = interpolate) +
    #geom_raster(interpolate = interpolate) +
    geom_vline( xintercept = 0,
                linetype = "dashed",
                linewidth = 2) +
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
#' @param chan_info channel info structure
#' @param sig_names names of signals in the data
#' @return Vector of channel positions
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

arrange_chans <- function(chan_info, sig_names) {

  # Use the sterographic projections to check which electrodes are midline
  midline <- is.true(chan_info$x == 0)

  # get labels of channels that are on the midline
  midline_labels <- chan_info$electrode[midline]
  midline_dist <- chan_info$y[midline]

  # Order them from back to front
  midline_labels <- midline_labels[sort(midline_dist,
                                        decreasing = TRUE,
                                        index.return = TRUE)$ix]

  # Pick out electrodes from the left hemisphere
  left <- is.true(chan_info$x < 0)
  left_labels <- chan_info$electrode[left]
  left_dist <- chan_info$y[left]
  left_labels <- left_labels[sort(left_dist,
                                  index.return = TRUE)$ix]

  # Pick out electrodes from the right hemisphere
  right <- is.true(chan_info$x > 0)
  right_labels <- chan_info$electrode[right]
  right_dist <- chan_info$y[right]
  right_labels <- right_labels[sort(right_dist,
                                    index.return = TRUE)$ix]

  all_labels <- c(left_labels,
                  midline_labels,
                  right_labels)

  elecs_found <- toupper(sig_names) %in% toupper(all_labels)

  new_ord <- match(toupper(c(all_labels,
                             sig_names[!elecs_found])),
                   toupper(sig_names))
  new_ord
}
