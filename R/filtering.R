#' Butterworth IIR filter
#'
#' Construct a Butterworth IIR filter and filter input data. This uses \code{signal::filt_filt}, which filters the signal twice to - once forwards, then again backwards)
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to be filtered.
#' @param low_freq Low passband edge.
#' @param high_freq High passband edge.
#' @param filter_order Order of the Butterworth filter.
#' @param srate Sampling rate of the signal.
#' @param plot_filt Plot filter characteristics
#' @param silent Turns off filtering messages.
#' @importFrom signal butter filtfilt freqz
#' @importFrom purrr map_df
#' @export
#'

iir_filt <- function(data, low_freq = NULL, high_freq = NULL, filter_order = 4,
                     srate = NULL, plot_filt = FALSE, silent = FALSE) {

  if (is.eeg_data(data)) {
    srate <- data$srate
    tmp_data <- data
    data <- data$signals
    eeg_dat <- TRUE
  } else if (is.null(srate)) {
    stop("sampling rate must be supplied.")
  } else {
    eeg_dat <- FALSE
  }

  if (any(filter_order < 2 | filter_order > 12)) {
    stop("Filter order should be between 2 and 12.")
  }

  if (is.null(low_freq)) {
    if (is.null(high_freq)) {
      stop("At least one frequency must be specified.")
    } else {
      filt_type <- "low"
      message(sprintf("Low-pass IIR filter at %.4g Hz", high_freq))
      W <- high_freq / (srate / 2)
    }
  } else if (is.null(high_freq)) {
    filt_type <- "high"
    message(sprintf("High-pass IIR filter at %.4g Hz", high_freq))
    W <- low_freq / (srate / 2)
    if (length(dim(data)) > 1) {
      data <- sweep(data, 2, colMeans(data))
    } else {
      data <- data - mean(data)
    }

  } else if (low_freq > high_freq) {
    filt_type <- "stop"
    message(sprintf("Band-stop IIR filter from %.4g-%.4g Hz",
                    low_freq, high_freq))
    W <- c(low_freq / (srate / 2), high_freq / (srate / 2))
    if (length(dim(data)) > 1) {
      data <- sweep(data, 2, colMeans(data))
    } else {
      data <- data - mean(data)
    }
  } else if (low_freq < high_freq) {
    filt_type <- "pass"
    message(sprintf("Band-pass IIR filter from %.4g-%.4g Hz",
                    low_freq, high_freq))
    W <- c(low_freq / (srate / 2), high_freq / (srate / 2))
    if (length(dim(data)) > 1) {
      data <- sweep(data, 2, colMeans(data))
    } else {
      data <- data - mean(data)
    }
  }

  #filtfilt filters twice, so effectively doubles filter_order
  filter_order <- round(filter_order / 2)
  filt_coef <- signal::butter(filter_order, W, type = filt_type)

  if (plot_filt) {
    signal::freqz(filt_coef, Fs = srate)
  }

  if (length(dim(data)) > 1) {
    if (eeg_dat) {
      tmp_data$signals <- purrr::map_df(as.list(data),
                                        ~signal::filtfilt(filt_coef, .))
      data <- tmp_data
    } else {
      data <- purrr::map_df(as.list(data), ~signal::filtfilt(filt_coef, .))
    }
  } else {
    data <- signal::filtfilt(filt_coef, data)
  }
  return(data)
}
