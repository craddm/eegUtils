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
#' Plot TFR objects
#'
#' @param data object of class eeg_tfr
#' @param electrode Electrode to plot
#' @param interpolate interpolation of raster
#' @param time_lim Time limits of plot
#' @noRd

plot_tfr <- function(data,
                     electrode,
                     interpolate = FALSE,
                     time_lim = NULL,
                     ...) {

  if (!class(data) == "eeg_tfr") {
    stop("Object of class eeg_tfr required.")
  }

  if (length(data$dimensions) == 4) {
    data$signals <- apply(data$signals,
                          c(1, 2, 3),
                          mean)
  }

  if (electrode %in% dimnames(data$signals)[[2]]) {
    aa <- as.data.frame.table(data$signals[, electrode, ],
                              stringsAsFactors = FALSE)
    names(aa) <- c("time", "frequency", "power")
    aa$time <- as.numeric(aa$time)
    aa$frequency <- as.numeric(aa$frequency)
    aa$power <- as.numeric(aa$power)
    #aa <- tidyr::gather(aa, frequency, power, -time)
    ggplot2::ggplot(aa, aes(x = time,
                            y = frequency,
                            fill = power)) +
      geom_raster(interpolate = interpolate) +
      labs(y = "Frequency (Hz)",
           x = "Time (s)",
           fill = "Power (a.u.)")
  } else {
    stop("Electrode not found.")
  }

}
