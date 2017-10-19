#' Construct a Butterworth IIR filter.
#'
#' @param data Data to be filtered.
#' @param low_freq Low passband edge.
#' @param high_freq High passband edge.
#' @param filter_order Order of the Butterworth filter.
#' @param srate Sampling rate of the signal.
#' @importFrom signal butter filtfilt
#'

iir_filt <- function(data, low_freq = NULL, high_freq = NULL, filter_order = 4, srate = NULL) {

  if (is.null(srate)) {
    error("sampling rate must be supplied.")
  }

  if (any(filter_order < 2 | filter_order > 10)) {
    error("Filter order should be between 2 and 10.")
  }

  if (is.null(low_freq)) {
    if (is.null(high_freq)) {
      error('At least one frequency must be specified.')
    } else {
      filt_type <- "low"
      message(sprintf('Low-pass IIR filter at %.4g Hz', low_freq))
      low_freq <- low_freq / (srate/2)
    }
  } else if (is.null(high_freq)) {
    filt_type <- "high"
    message(sprintf('High-pass IIR filter at %.4g Hz', high_freq))
    high_freq <- high_freq / (srate/2)
    data <- data - mean(data)
  } else if (low_freq > high_freq) {
    filt_type <- "stop"
    message(sprintf('Band-stop IIR filter from %.4g-%.4g Hz', low_freq, high_freq))
    low_freq <- low_freq / (srate/2)
    high_freq <- high_freq / (srate/2)
    data <- data - mean(data)
  } else if (low_freq < high_freq) {
    filt_type <- "pass"
    message(sprintf('Band-pass IIR filter from %.4g-%.4g Hz', low_freq, high_freq))
    low_freq <- low_freq / (srate/2)
    high_freq <- high_freq / (srate/2)
    data <- data - mean(data)
  }

  filter_order <- filter_order/2 #filtfilt filters twice, so effectively doubles filter_order
  filt_coef <- signal::butter(filter_order, c(low_freq, high_freq), type = filt_type)
  signal::filtfilt(filt_coef, data)
}

