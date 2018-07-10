#' Plot power spectral density
#'
#' @param data Object of class \code{eeg_epochs}
#' @param ... Additional parameters
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export
plot_psd <- function(data, freq_range = NULL, ...) {
  UseMethod("plot_psd", data)
}

#' @param freq_range Vector of lower and upper frequencies to plot. (e.g. c(1, 40))
#' @describeIn plot_psd Plot PSDs for \code{eeg_epochs}.
#' @export
plot_psd.eeg_epochs <- function(data, freq_range = NULL, ...) {

  psd_out <- compute_psd(data, keep_trials = FALSE, n_fft = 512)

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
    theme_bw()
}

#' @describeIn plot_psd Plot PSDs for \code{eeg_epochs}.
#' @export
plot_psd.eeg_data <- function(data, freq_range = NULL, ...) {

  psd_out <- compute_psd(data, n_fft = 512)

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
    theme_bw()
}

#' plot tfr objects
#' @param data object of class eeg_tfr
#' @noRd
plot_tfr <- function(data, electrode, ...) {
  if (!class(data) == "eeg_tfr") {
    stop("Object of class eeg_tfr required.")
  }
  #length(dim(data$signals))

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
      geom_raster(interpolate = TRUE)
  } else {
    stop("Electrode not found.")
  }

}
