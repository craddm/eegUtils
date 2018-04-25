#' Butterworth IIR filter
#'
#' Construct a Butterworth IIR filter and filter input data. This uses
#' \code{signal::filt_filt}, which filters the signal twice to - once forwards,
#' then again backwards).
#'
#' low_freq and high_freq are passband edges. Pass low freq or high freq alone
#' to perform high-pass or low-pass filtering respectively. For band-pass or
#' band-stop filters, pass both low_freq and high_freq.
#'
#' If low_freq < high_freq, bandpass filtering is performed.
#'
#' If low_freq > high_freq, bandstop filtering is performed.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to be filtered.
#' @param ... Parameters passed to S3 methods
#' @export

iir_filt <- function(data, ...) {
  UseMethod("iir_filt", data)
}

#' @export
iir_filt.default <- function(data, ...) {
  stop("Function works on data.frames, eeg_data, and eeg_epochs.")
}

#' @param low_freq Low passband edge.
#' @param high_freq High passband edge.
#' @param filter_order Order of the Butterworth filter.
#' @param srate Sampling rate of the signal.
#' @param silent Turns off filtering messages.
#' @importFrom signal butter filtfilt freqz
#' @importFrom purrr map_df
#' @export
#' @describeIn iir_filt Filter a data frame
#'

iir_filt.data.frame <- function(data, low_freq = NULL, high_freq = NULL,
                                filter_order = 4, srate,
                                silent = FALSE, ...) {
  if (missing(srate)) {
    stop("sampling rate must be supplied.")
  }
  data <- run_iir(data, low_freq, high_freq, filter_order, srate, silent = silent)
  data
}

#' @export
#' @describeIn iir_filt Filter eeg_data

iir_filt.eeg_data <- function(data, low_freq = NULL, high_freq = NULL,
                              filter_order = 4, silent = FALSE, ...) {

  if (!is.null(data$reference)) {
    data$signals["ref_data"] <- data$reference$ref_data
  }

  data$signals <- run_iir(data$signals, low_freq, high_freq, filter_order,
                          data$srate, silent = silent)

  if (!is.null(data$reference)) {
    data$reference$ref_data <- data$signals["ref_data"]
    data$signals["ref_data"] <- NULL
  }

  return(data)
}

#' @export
#' @describeIn iir_filt Filter eeg_epochs
iir_filt.eeg_epochs <- function(data, low_freq = NULL, high_freq = NULL,
                                filter_order = 4,
                                silent = FALSE, ...) {

  if (!is.null(data$reference)) {
    data$signals["ref_data"] <- data$reference$ref_data
  }

  data$signals$epoch <- data$timings$epoch
  data$signals <- run_iir(data$signals, low_freq, high_freq, filter_order,
                          data$srate, silent = silent)
  data$signals["epoch"] <- NULL

  if (!is.null(data$reference)) {
    data$reference$ref_data <- data$signals["ref_data"]
    data$signals["ref_data"] <- NULL
  }

  data
}


#' Internal function for running IIR filtering
#'
#' @param data Data to be filtered
#' @param low_freq Low passband edge.
#' @param high_freq High passband edge.
#' @param filter_order Order of the Butterworth filter.
#' @param srate Sampling rate of the signal.
#' @param silent Turns off filtering messages.
#' @importFrom dplyr group_by
#' @importFrom purrr map_df
#' @importFrom signal filtfilt butter

run_iir <- function(data, low_freq = NULL, high_freq = NULL, filter_order,
                    srate, silent = FALSE) {

  if (filter_order < 2 || filter_order > 20) {
    stop("Filter order should be between 2 and 20.")
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
    message(sprintf("High-pass IIR filter at %.4g Hz", low_freq))
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

  #filtfilt filters twice, so effectively doubles filter_order - we half it here
  #so that it corresponds to the expectation of the user
  filter_order <- round(filter_order / 2)
  filt_coef <- signal::butter(filter_order, W, type = filt_type)

  if ("epoch" %in% names(data)) {
    data <- dplyr::group_by(data, epoch)
  }

  if ("electrode" %in% names(data)) {
    data <- dplyr::group_by(data, electrode, add = TRUE)
  }

  if (length(dim(data)) > 1) {
    data <- purrr::map_df(as.list(data), ~signal::filtfilt(filt_coef, .))
    } else {
    data <- signal::filtfilt(filt_coef, data)
  }
  data
}
