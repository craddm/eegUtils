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
#' Note that the signal is first zero-meaned using either channel means or
#' by-channel epoch means.
#' @examples
#' plot_psd(iir_filt(demo_epochs, low_freq = 1, high_freq = 30))
#' plot_psd(iir_filt(demo_epochs, low_freq = 12, high_freq = 8))
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

iir_filt.data.frame <- function(data,
                                low_freq = NULL,
                                high_freq = NULL,
                                filter_order = 4,
                                srate,
                                silent = FALSE,
                                ...) {
  if (missing(srate)) {
    stop("sampling rate must be supplied.")
  }
  data <- run_iir(data,
                  low_freq,
                  high_freq,
                  filter_order,
                  srate,
                  silent = silent)
  data
}

#' @export
#' @describeIn iir_filt Filter \code{eeg_data} objects

iir_filt.eeg_data <- function(data,
                              low_freq = NULL,
                              high_freq = NULL,
                              filter_order = 4,
                              silent = FALSE,
                              ...) {

  data$signals <- run_iir(data$signals,
                          low_freq,
                          high_freq,
                          filter_order,
                          data$srate,
                          silent = silent)
  data
}

#' @export
#' @describeIn iir_filt Filter \code{eeg_epochs} objects.
iir_filt.eeg_epochs <- function(data,
                                low_freq = NULL,
                                high_freq = NULL,
                                filter_order = 4,
                                silent = FALSE,
                                ...) {

  data$signals$epoch <- data$timings$epoch
  data$signals <- run_iir(data$signals,
                          low_freq,
                          high_freq,
                          filter_order,
                          data$srate,
                          silent = silent)
  data$signals["epoch"] <- NULL
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
#' @keywords internal

run_iir <- function(data,
                    low_freq = NULL,
                    high_freq = NULL,
                    filter_order,
                    srate,
                    silent = FALSE) {

  if (filter_order < 2 || filter_order > 20) {
    stop("Filter order should be between 2 and 20.")
  }

  if (length(low_freq) > 1 | length(high_freq) > 1) {
    stop("Only one number should be passed to low_freq or high_freq")
  }

  if (is.null(low_freq) & is.null(high_freq)) {
    stop("At least one frequency must be specified.")
  }

  if (is.null(low_freq)) {
      filt_type <- "low"
      message(sprintf("Low-pass IIR filter at %.4g Hz", high_freq))
      W <- high_freq / (srate / 2)
    } else if (is.null(high_freq)) {
    filt_type <- "high"
    message("High-pass IIR filter at ", low_freq," Hz")
    W <- low_freq / (srate / 2)

    if (length(dim(data)) > 1) {
      data <- sweep(data, 2, colMeans(data))
    } else {
      data <- data - mean(data)
    }

  } else if (low_freq > high_freq) {
    filt_type <- "stop"
    message(sprintf("Band-stop IIR filter from %.4g-%.4g Hz",
                    high_freq, low_freq))
    W <- c(high_freq / (srate / 2), low_freq / (srate / 2))

    if (length(dim(data)) > 1) {
      data <- sweep(data, 2, colMeans(data))
    } else {
      data <- data - mean(data)
    }
  } else if (low_freq < high_freq) {
    filt_type <- "pass"
    message(sprintf("Band-pass IIR filter from %.4g-%.4g Hz",
                    low_freq, high_freq))
    W <- c(low_freq / (srate / 2),
           high_freq / (srate / 2))

    if (length(dim(data)) > 1) {
      data <- sweep(data, 2, colMeans(data))
    } else {
      data <- data - mean(data)
    }
  }

  #filtfilt filters twice, so effectively doubles filter_order - we half it here
  #so that it corresponds to the expectation of the user
  filter_order <- round(filter_order / 2)
  filt_coef <- signal::butter(filter_order,
                              W,
                              type = filt_type)

  data <- data.table(data)

  if ("epoch" %in% names(data)) {
    data <- data[, lapply(.SD, function(x) signal::filtfilt(filt_coef, x)), by = epoch]
  } else {
    data <- data[, lapply(.SD, function(x) signal::filtfilt(filt_coef, x))]
  }

  # if (length(dim(data)) > 1) {
  #   data <- furrr::future_map_dfr(as.list(data),
  #                         ~signal::filtfilt(filt_coef, .))
  #   } else {
  #   data <- signal::filtfilt(filt_coef,
  #                            data)
  # }
  tibble::as_tibble(data)
}


#' Gaussian filter
#'
#' Gaussian filtering in the frequency domain.
#'
#' @param data Data to be filtered
#' @param srate Sampling rate of the data
#' @param freq Peak frequency of the filter
#' @param fwhm Standard deviation of the filter
#' @noRd

gauss_filter <- function(data,
                         srate,
                         freq,
                         fwhm) {
  hz <- seq(0, srate,
            length.out = nrow(data))
  s <- fwhm * (2 * pi - 1) / (4 * pi)
  x <- hz - freq
  fx <- exp(-.5 * (x / s) ^2)
  fx <- fx / max(fx)
  filt_sig <- apply(data, 2, function(x) {
    2 * Re(fft(fft(x) / srate * fx,
               inverse = TRUE))})
  as.data.frame(filt_sig)
}
