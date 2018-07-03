#' @param data TFR transformed data
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param freqs vector of frequencies
#' @param dimensions List of which dimension is which
#' @noRd
eeg_tfr <- function(data,
                    srate,
                    events,
                    chan_info = NULL,
                    timings = NULL,
                    freqs,
                    dimensions) {

  value <- list(signals = data,
                srate = srate,
                events = events,
                chan_info = chan_info,
                timings = timings,
                freqs = freqs,
                dimensions = dimensions)
  class(value) <- "eeg_tfr"
  value
}


#' Function to create an object of class \code{eeg_psd}
#'
#' @param data PSD transformed data
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param freqs vector of frequencies
#' @param dimensions List of which dimension is which
#' @noRd
eeg_psd <- function(data,
                    srate,
                    events,
                    chan_info = NULL,
                    timings = NULL,
                    freqs,
                    dimensions) {

  value <- list(signals = data,
                srate = srate,
                chan_info = chan_info,
                timings = timings,
                freqs = freqs
                )
  class(value) <- "eeg_tfr"
  value
}
