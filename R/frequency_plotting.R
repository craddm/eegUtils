#' Plot Power Spectral Density
#'
#' Calculate and plot the PSD for \code{eeg_*} objects. Output units are dB. The
#' PSD is calculated using Welch's method.
#'
#' Welch's method splits the data into multiple segments and then averages over
#' those segments. For epoched data, Welch's FFT is calculated separately for
#' each trial.
#'
#' Specific parameters such as the number of FFT points and the amount of
#' overlap between segments can be passed to Welch's FFT
#'
#' @examples
#'  plot_psd(demo_epochs)
#'  plot_psd(demo_epochs, seg_length = 256)
#' @param data Object of class \code{eeg_epochs}, \code{eeg_data}, or \code{eeg_ICA}.
#' @param freq_range Vector of lower and upper frequencies to plot. (e.g. c(1,
#'   40))
#' @param ... Additional parameters.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export

plot_psd <- function(data, freq_range = NULL, ...) {
  UseMethod("plot_psd", data)
}

#' @param n_fft Number of points to use for the underlying FFTs. Defaults to 512
#'   for \code{eeg_epochs} or minimum of 2048 or the signal length for
#'   \code{eeg_data}.
#' @param noverlap Amount of overlap between segments, in sampling points.
#'   Defaults to 50\%.
#' @param seg_length Length of individual segments. Defaults to n_fft. Must be <= n_fft.
#' @describeIn plot_psd Plot PSD for \code{eeg_epochs}.
#' @export
plot_psd.eeg_epochs <- function(data,
                                freq_range = NULL,
                                n_fft = 512,
                                seg_length = NULL,
                                noverlap = NULL,
                                ...) {

  psd_out <- compute_psd(data,
                         keep_trials = FALSE,
                         n_fft = n_fft,
                         seg_length = seg_length,
                         noverlap = noverlap)

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

#' @describeIn plot_psd Plot PSD for \code{eeg_data}.
#' @export
plot_psd.eeg_data <- function(data,
                              freq_range = NULL,
                              n_fft = 2048,
                              noverlap = NULL,
                              seg_length = NULL,
                              ...) {

  psd_out <- compute_psd(data,
                         n_fft = n_fft,
                         seg_length = NULL,
                         noverlap = NULL)

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

#' @param components Which components to compute the PSD for. Defaults to all.
#' @describeIn plot_psd Plot PSD for \code{eeg_ICA} objects
#' @export
plot_psd.eeg_ICA <- function(data,
                             freq_range = NULL,
                             components = NULL,
                             seg_length = NULL,
                             noverlap = NULL,
                             n_fft = 512,
                             ...) {

  if (!is.null(components)) {
    data <- select_elecs(data,
                         components)
  }

  psd_out <- compute_psd(data,
                         n_fft = n_fft,
                         seg_length = seg_length,
                         noverlap = noverlap,
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

#' @describeIn plot_psd Plot PSD for \code{data.frame}s.
#' @importFrom dplyr select group_by summarise_all
#' @export
plot_psd.data.frame <- function(data,
                                freq_range = NULL,
                                ...) {

  if ("epoch" %in% names(data)) {
    data <- dplyr::select(data, -epoch)
    data <- dplyr::group_by(data, frequency)
    data <- dplyr::summarise_all(data, mean)
  }
  data <- tidyr::gather(data,
                        electrode,
                        power,
                        -frequency)

  data$power <- 10 * log10(data$power)
  ggplot(data,
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
#' Creates a time-frequency plot of an \code{eeg_tfr} object. The plot has time
#' on the x-axis and frequency on the y-axis. If no electrode is supplied, it
#' will average over all electrodes.
#'
#' Various different baseline options can be applied here (e.g. "db" for
#' decibels, "pc" for percent change, "divide" for division; see
#' \code{rm_baseline} for details).
#'
#' @param data Object of class \code{eeg_tfr}
#' @param electrode Electrode to plot. If none is supplied, averages over all
#'   electrodes.
#' @param time_lim Time limits of plot.
#' @param freq_range Vector of two numbers. (e.g. c(8, 40)).
#' @param baseline Baseline period
#' @param baseline_type baseline correction to apply. Defaults to "none".
#' @param fill_lims Custom colour scale (i.e. range of power). e.g. c(-5, 5).
#' @param interpolate Interpolation of raster for smoother plotting.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @importFrom purrr partial
#' @import ggplot2
#' @seealso [rm_baseline()]
#' @export

plot_tfr <- function(data,
                     electrode = NULL,
                     time_lim = NULL,
                     freq_range = NULL,
                     baseline_type = "none",
                     baseline = NULL,
                     fill_lims = NULL,
                     interpolate = FALSE) {

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
    data$dimensions <- data$dimensions[1:3]
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
    fill_lims <- c(min(data$signals, na.rm = TRUE),
                   max(data$signals, na.rm = TRUE))
    if (all(fill_lims > 0)) {
      fill_lims <- c(0, max(abs(fill_lims)))
    } else {
      fill_lims <- max(abs(fill_lims))
      fill_lims <- c(-fill_lims, fill_lims)
    }
  }

  # Use purrr::partial to save copy pasting the whole thing in every
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
