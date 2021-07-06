#' Plot Power Spectral Density
#'
#' Calculate and plot the PSD for `eeg_*` objects. Output units are dB. The
#' PSD is calculated using Welch's method.
#'
#' Welch's method splits the data into multiple segments and then averages over
#' those segments. For epoched data, Welch's FFT is calculated separately for
#' each trial.
#'
#' Specific parameters such as the number of FFT points and the amount of
#' overlap between segments can be passed to Welch's FFT.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @examples
#'  plot_psd(demo_epochs)
#'  plot_psd(demo_epochs, seg_length = 256)
#' @param data Object of class `eeg_epochs`, `eeg_data`, or `eeg_ICA`.
#' @param freq_range Vector of lower and upper frequencies to plot. (e.g. c(1,
#'   40))
#' @param ... Additional parameters.
#' @return A ggplot object.
#' @export

plot_psd <- function(data, freq_range = NULL, ...) {
  UseMethod("plot_psd", data)
}

#' @param n_fft Number of points to use for the underlying FFTs. Defaults to 256
#'   for `eeg_epochs` or minimum of 2048 or the signal length for `eeg_data`.
#' @param noverlap Amount of overlap between segments, in sampling points.
#'   Defaults to 50%.
#' @param seg_length Length of individual segments. Defaults to n_fft. Must be
#'   <= n_fft.
#' @param demean Remove epoch means before FFT.
#' @param keep_trials Whether to keep trial information in the output or average
#'   over all trials
#' @describeIn plot_psd Plot PSD for `eeg_epochs`.
#' @export
plot_psd.eeg_epochs <- function(data,
                                freq_range = NULL,
                                n_fft = 256,
                                seg_length = NULL,
                                noverlap = NULL,
                                demean = TRUE,
                                keep_trials = TRUE,
                                ...) {

  psd_out <- compute_psd(data,
                         keep_trials = keep_trials,
                         n_fft = n_fft,
                         seg_length = seg_length,
                         noverlap = noverlap,
                         demean = demean)

  create_psd_plot(psd_out,
                  freq_range,
                  chan_names = channel_names(data))
}

#' @describeIn plot_psd Plot PSD for `eeg_data`.
#' @export
plot_psd.eeg_data <- function(data,
                              freq_range = NULL,
                              n_fft = 2048,
                              noverlap = NULL,
                              seg_length = NULL,
                              ...) {

  psd_out <- compute_psd(data,
                         n_fft = n_fft,
                         seg_length = seg_length,
                         noverlap = noverlap)

  create_psd_plot(psd_out,
                  freq_range,
                  chan_names = channel_names(data))

}

#' @param components Which components to compute the PSD for. Defaults to all.
#' @describeIn plot_psd Plot PSD for `eeg_ICA` objects
#' @export
plot_psd.eeg_ICA <- function(data,
                             freq_range = NULL,
                             components = NULL,
                             seg_length = NULL,
                             noverlap = NULL,
                             n_fft = 256,
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

  create_psd_plot(psd_out,
                  freq_range,
                  chan_names = channel_names(data))
}

#' @describeIn plot_psd Plot PSD for `data.frame`s.
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

#' @describeIn plot_psd Plot PSD for `eeg_evoked` objects
#' @export
plot_psd.eeg_evoked <- function(data,
                                freq_range = NULL,
                                n_fft = 256,
                                seg_length = NULL,
                                noverlap = NULL,
                                keep_trials = TRUE,
                                ...) {

  psd_out <- compute_psd(data,
                         keep_trials = keep_trials,
                         n_fft = n_fft,
                         seg_length = seg_length,
                         noverlap = noverlap)

  create_psd_plot(psd_out,
                  freq_range,
                  chan_names = channel_names(data))
}

#' @describeIn plot_psd Plot PSD for `eeg_group` objects is not currently supported
#' @export
plot_psd.eeg_group <- function(data,
                               freq_range = NULL,
                               n_fft = 256,
                               seg_length = NULL,
                               noverlap = NULL,
                               demean = TRUE,
                               ...) {
  stop("Cannot currently plot_psd for eeg_group objects.")

}


#' Create a PSD plot
#'
#' @param psd_out PSD to plot.
#' @param freq_range Frequency range to plot.
#' @return ggplot showing power spectral density.
#' @keywords internal
create_psd_plot <- function(psd_out,
                            freq_range,
                            chan_names) {

  if (!is.null(freq_range)) {
    if (length(freq_range) < 2 | length(freq_range) > 2) {
      message("freq_range must be a vector of length 2. Displaying all frequencies.")
    } else {
      rows <- psd_out$frequency >= freq_range[[1]] &
        psd_out$frequency <= freq_range[[2]]
      psd_out <- psd_out[rows, ]
    }
  }

  psd_out <-
    tidyr::pivot_longer(psd_out,
                        cols = chan_names,
                        names_to = "electrode",
                        values_to = "power")
  ggplot(psd_out,
         aes(x = frequency,
             y = 10 * log10(power),
             colour = electrode)) +
    stat_summary(geom = "line",
                 fun = mean) +# geom_line() +
    theme_bw() +
    ylab(expression(paste(mu, V^2, "/ Hz(dB)"))) +
    xlab("Frequency (Hz)") +
    scale_x_continuous(expand = c(0, 0))
}

#' Time-frequency plot
#'
#' Creates a time-frequency plot of an `eeg_tfr` object. The plot has time
#' on the x-axis and frequency on the y-axis. If no electrode is supplied, it
#' will average over all electrodes.
#'
#' Various different baseline options can be applied here (e.g. "db" for
#' decibels, "pc" for percent change, "divide" for division; see
#' `rm_baseline` for details).
#'
#' @param data Object of class `eeg_tfr`
#' @param electrode Electrode to plot. If none is supplied, averages over all
#'   electrodes.
#' @param time_lim Time limits of plot.
#' @param freq_range Vector of two numbers. (e.g. c(8, 40)).
#' @param baseline Baseline period
#' @param baseline_type baseline correction to apply. Defaults to "none".
#' @param fill_lims Custom colour scale (i.e. range of power). e.g. c(-5, 5).
#' @param interpolate Interpolation of raster for smoother plotting.
#' @param na.rm Remove NA values silently (TRUE) or with a warning (FALSE).
#'   Defaults to TRUE.
#' @param fun.data Statistical function to use for averaging over
#'   electrodes/conditions. Defaults to `mean`.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @importFrom purrr partial
#' @import ggplot2
#' @return A `ggplot`
#' @seealso [rm_baseline()]
#' @export

plot_tfr <- function(data,
                     electrode = NULL,
                     time_lim = NULL,
                     freq_range = NULL,
                     baseline_type = "none",
                     baseline = NULL,
                     fill_lims = NULL,
                     interpolate = FALSE,
                     na.rm = TRUE,
                     fun.data = mean) {

  if (!inherits(data,"eeg_tfr")) {
    stop("Object of class eeg_tfr required.")
  }

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  if (!is.null(electrode)) {
    data <- select_elecs(data, electrode)
  }

  if (!is.null(freq_range)) {
    data_freqs <- as.numeric(dimnames(data$signals)[["frequency"]])
    data_freqs <- (data_freqs >= freq_range[1] & data_freqs <= freq_range[2])
    data$signals <- data$signals[, , , data_freqs, drop = FALSE]
  }

  if (identical(data$freq_info$output, "fourier")) {
    message("Data is complex Fourier coefficients, converting to power.")
    data <- convert_tfr(data,
                        "power")
    data$freq_info$output <- "power"
  }

  if (!inherits(data, c("tfr_average",
                        "eeg_group"))) {
    data <- eeg_average(data)
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
    if (identical(data$freq_info$baseline, "none")) {
      fill_lims <- c(0, NA)
    } else{
      fill_lims <- c(NA, NA)
    }
  }

  # Use purrr::partial to save copy pasting the whole thing in every
  fill_dist <- purrr::partial(scale_fill_distiller,
                              palette = "RdBu",
                              limits = fill_lims,
                              oob = scales::squish)

  fill_colour <-
    switch(data$freq_info$baseline,
           "none" = ggplot2::scale_fill_viridis_c(limits = fill_lims,
                                                  oob = scales::squish),
           "absolute" = fill_dist(),
           "db" = fill_dist(),
           "divide" = fill_dist(),
           "ratio" = fill_dist(),
           "pc" = fill_dist(),
           ggplot2::scale_fill_viridis_c())

  data <- as.data.frame(data,
                        long = TRUE)

  tfr_plot <-
    ggplot2::ggplot(data,
                    aes(x = time,
                        y = frequency,
                        fill = power)) +
    stat_summary_by_fill(fun.data = mean,
                         na.rm = na.rm,
                         interpolate = interpolate) +
    labs(x = "Time (s)",
         y = "Frequency (Hz)",
         fill = fill_lab) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    fill_colour

  if (is.unsorted(diff(unique(data$frequency)))) {
    tfr_plot <-
      tfr_plot +
      scale_y_continuous(expand = c(0, 0))
  } else {
    tfr_plot <-
      tfr_plot +
      scale_y_log10(expand = c(0, 0))
  }



  tfr_plot
}
