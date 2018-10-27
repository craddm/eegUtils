#' Plot power spectral density
#'
#' Calculate and plot the PSD for \code{eeg_*} objects. Output units are dB. 512
#' points are used for FFTs.
#'
#' @param data Object of class \code{eeg_epochs}
#' @param freq_range Vector of lower and upper frequencies to plot. (e.g. c(1,
#'   40))
#' @param ... Additional parameters
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export
plot_psd <- function(data, freq_range = NULL, ...) {
  UseMethod("plot_psd", data)
}

#' @describeIn plot_psd Plot PSDs for \code{eeg_epochs}.
#' @export
plot_psd.eeg_epochs <- function(data,
                                freq_range = NULL,
                                ...) {

  psd_out <- compute_psd(data,
                         keep_trials = FALSE,
                         n_fft = 512)

  if (!is.null(freq_range)) {
    if (length(freq_range) < 2 | length(freq_range) > 2) {
      message("freq_range must be a vector of length 2. Displaying all frequencies.")
    } else {
      rows <- psd_out$frequency >= freq_range[[1]] &
        psd_out$frequency <= freq_range[[2]]
      psd_out <- psd_out[rows, ]
    }
  }

  psd_out <- tidyr::gather(psd_out,
                           electrode,
                           power,
                           -frequency)
  psd_out$power <- 10 * log10(psd_out$power)
  ggplot(psd_out,
         aes(x = frequency,
             y = power,
             colour = electrode)) +
    geom_line() +
    theme_bw() +
    ylab("Decibels (10 * log10(uV^2 / Hz)") +
    xlab("Frequency (Hz)")
}

#' @describeIn plot_psd Plot PSDs for \code{eeg_epochs}.
#' @export
plot_psd.eeg_data <- function(data,
                              freq_range = NULL,
                              ...) {

  psd_out <- compute_psd(data,
                         n_fft = 512)

  if (!is.null(freq_range)) {
    if (length(freq_range) < 2 | length(freq_range) > 2) {
      message("freq_range must be a vector of length 2. Displaying all frequencies.")
    } else {
      rows <- psd_out$frequency >= freq_range[[1]] &
        psd_out$frequency <= freq_range[[2]]
      psd_out <- psd_out[rows, ]
    }
  }

  psd_out <- tidyr::gather(psd_out,
                           electrode,
                           power,
                           -frequency)
  psd_out$power <- 10 * log10(psd_out$power)
  ggplot(psd_out,
         aes(x = frequency,
             y = power,
             colour = electrode)) +
    geom_line() +
    theme_bw() +
    ylab("Decibels (10 * log10(uV^2 / Hz)") +
    xlab("Frequency (Hz)")
}

#'@param components Which components to compute the PSD for
#'@describeIn plot_psd Plot PSD for \code{eeg_ICA} objects
#'@export
plot_psd.eeg_ICA <- function(data,
                             freq_range = NULL,
                             components = NULL,
                             ...) {

  if (!is.null(components)) {
    data <- select_elecs(data,
                         components)
  }

  psd_out <- compute_psd(data,
                         n_fft = 512,
                         keep_trials = FALSE)

  if (!is.null(freq_range)) {
    if (length(freq_range) < 2 | length(freq_range) > 2) {
      message("freq_range must be a vector of length 2. Displaying all frequencies.")
    } else {
      rows <- psd_out$frequency >= freq_range[[1]] &
        psd_out$frequency <= freq_range[[2]]
      psd_out <- psd_out[rows, ]
    }
  }

  psd_out <- tidyr::gather(psd_out,
                           electrode,
                           power,
                           -frequency)
  psd_out$power <- 10 * log10(psd_out$power)
  ggplot(psd_out,
         aes(x = frequency,
             y = power,
             colour = electrode)) +
    geom_line() +
    theme_bw() +
    ylab("Decibels (10 * log10(uV^2 / Hz)") +
    xlab("Frequency (Hz)")

}
#' Time-frequency plot
#'
#' Create time-frequency plot of an \code{eeg_tfr} object.
#'
#' Various different baseline options can be applied.
#'
#' @param data Object of class \code{eeg_tfr}
#' @param electrode Electrode to plot. If none is supplied, averages over all electrodes.
#' @param time_lim Time limits of plot.
#' @param freq_range Vector of two numbers. (e.g. c(8, 40)).
#' @param baseline Baseline period
#' @param baseline_type baseline correction to apply. Defaults to "none".
#' @param fill_lims Custom colour scale (i.e. range of power). e.g. c(-5, 5).
#' @param interpolate Interpolation of raster for smoother plotting.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @importFrom purrr partial
#' @import ggplot2
#' @export

plot_tfr <- function(data,
                     electrode = NULL,
                     time_lim = NULL,
                     freq_range = NULL,
                     baseline_type = "none",
                     baseline = NULL,
                     fill_lims = NULL,
                     interpolate = FALSE,
                     ...) {

  if (!class(data) == "eeg_tfr") {
    stop("Object of class eeg_tfr required.")
  }

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  if (!is.null(electrode)) {
    data <- select_elecs(data, electrode)
  }

  if (!is.null(freq_range)) {
    data_freqs <- as.numeric(dimnames(data$signals)[[3]])
    data_freqs <- (data_freqs >= freq_range[1] & data_freqs <= freq_range[2])
    data$signals <- data$signals[, , data_freqs, drop = FALSE]
  }

  if (length(data$dimensions) == 4) {
    data$signals <- apply(data$signals,
                          c(1, 2, 3),
                          mean)
  }

  if (baseline_type != "none") {
    data <- rm_baseline(data,
                        time_lim = baseline,
                        type = baseline_type)
  }

  fill_lab <-
    switch(data$freq_info$baseline,
           "none" = "Power (a.u.)",
           "db" = "Power (dB)",
           "divide" = "Relative power (%)",
           "ratio" = "Power ratio",
           "pc" = "Percent change (%)",
           "absolute" = "Power (a.u.)",
           "Power (a.u.)")

  if (is.null(fill_lims)) {
    fill_lims <- abs(c(min(data$signals, na.rm = TRUE),
                       max(data$signals, na.rm = TRUE)))
    fill_lims <- max(fill_lims)
    fill_lims <- c(-fill_lims, fill_lims)
  }

  fill_dist <- purrr::partial(scale_fill_distiller,
                              palette = "RdBu",
                              limits = fill_lims,
                              oob = scales::squish)

  fill_colour <-
    switch(data$freq_info$baseline,
           "none" = scale_fill_viridis_c(limits = fill_lims,
                                         oob = scales::squish),
           "absolute" = fill_dist(),
           "db" = fill_dist(),
           "divide" = fill_dist(),
           "ratio" = fill_dist(),
           "pc" = fill_dist(),
           scale_fill_viridis_c())

  data <- as.data.frame(data, long = TRUE)

  tfr_plot <-
    ggplot2::ggplot(data,
                    aes(x = time,
                        y = frequency,
                        fill = power)) +
    geom_raster(interpolate = interpolate) +
    labs(x = "Time (s)",
         y = "Frequency (Hz)",
         fill = fill_lab) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    fill_colour

  tfr_plot
}
