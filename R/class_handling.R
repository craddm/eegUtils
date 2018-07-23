#' Object creator for eeg_tfr objects.
#'
#' @param data TFR transformed data
#' @param srate Sampling rate in Hz.
#' @param events Event tables
#' @param chan_info Standard channel information.
#' @param reference Reference information
#' @param timings Timing information.
#' @param freq_info Frequencies and other useful information
#' @param dimensions List of which dimension is which
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @keywords internal
eeg_tfr <- function(data,
                    srate,
                    events,
                    chan_info = NULL,
                    reference,
                    timings = NULL,
                    freq_info,
                    dimensions) {

  value <- list(signals = data,
                srate = srate,
                events = events,
                chan_info = chan_info,
                reference = reference,
                timings = timings,
                freq_info = freq_info,
                dimensions = dimensions)
  class(value) <- "eeg_tfr"
  value
}

#' Check if object is of class eeg_tfr
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object to check.
#'
is.eeg_tfr <- function(x) inherits(x, "eeg_tfr")


#' Function to create an object of class eeg_psd
#'
#' @param data PSD transformed data
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param freqs vector of frequencies
#' @param dimensions List of which dimension is which
#' @author Matt Craddock \email{matt@@mattcraddock.com}
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
