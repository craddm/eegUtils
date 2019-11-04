#' Compute power spectral density
#'
#' \code{compute_psd} returns the PSD calculated using Welch's method for every
#' channel in the data. The output is in  microvolts ^2 / Hz. If the object has
#' multiple epochs, it will perform Welch's FFT separately for each epoch and
#' then average them afterwards.
#'
#' Welch's FFT splits the data into multiple segments, calculates the FFT
#' separately for each segment, and then averages over segments. Each segment is
#' windowed with a Hanning window to counter spectral leakage. For epoched data,
#' Welch's FFT is calculated separately for each trial.
#'
#' The number of sampling points used for the FFT can be specified using n_fft.
#' n_fft defaults to 256 sampling points for \code{eeg_epochs} data, or the
#' minimum of 2048 or the length of the signal for continuous \code{eeg_data}.
#'
#' \code{seg_length} defaults to be \code{n_fft}, and must be less than or equal
#' to it.
#'
#' \code{noverlap} specifies the amount of overlap between windows in sampling
#' points. If not specified, it defaults to 50\% overlap between segments.
#'
#' @examples
#' compute_psd(demo_epochs)
#' compute_psd(demo_epochs, n_fft = 256, seg_length = 128)
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to be plotted. Accepts objects of class \code{eeg_data}
#' @param ... any further parameters passed to specific methods
#' @return Currently, a data frame with the PSD for each channel separately.
#' @export

compute_psd <- function(data, ...) {
  UseMethod("compute_psd", data)
}

#' @param n_fft Length of FFT to be calculated in sampling points. See details.
#' @param seg_length Length of rolling data segments. Defaults to \code{n_fft}.
#'   Must be <= \code{n_fft}.
#' @param noverlap Number of (sampling) points of overlap between segments. Must
#'   be <= \code{seg_length}.
#' @param method Defaults to "Welch". No other method currently implemented.
#' @param demean Remove channel/epoch means. TRUE by default.
#' @param verbose Print informative messages. TRUE by default.
#' @describeIn compute_psd Compute PSD for an \code{eeg_data} object
#' @export

compute_psd.eeg_data <- function(data,
                                 seg_length = NULL,
                                 noverlap = NULL,
                                 n_fft = NULL,
                                 method = "Welch",
                                 demean = TRUE,
                                 verbose = TRUE,
                                 ...) {

  if (demean) {
    data <- rm_baseline(data,
                        verbose = verbose)
  }

  srate <- data$srate

  if (is.null(n_fft)) {
    n_fft <- min(2048, c(nrow(data$signals)))
  }

  if (is.null(seg_length)) {
    seg_length <- n_fft
  }

  if (seg_length > n_fft) {
    stop("seg_length cannot be greater than n_fft")
  }

  if (is.null(noverlap)) {
    noverlap <- seg_length %/% 2
  } else if (noverlap >= seg_length) {
    stop("noverlap should not be larger than seg_length.")
  }

  if (method == "Welch") {
    final_output <- welch_fft(data$signals,
                              seg_length,
                              noverlap = noverlap,
                              n_fft = n_fft,
                              srate = srate,
                              n_sig = nrow(data$signals))

  }  else {
    stop("Welch is the only available method at this time.")
  }
  final_output
}

#' @param keep_trials Include FFT for every trial in output, or average over them if FALSE.
#' @describeIn compute_psd Compute PSD for an \code{eeg_epochs} object
#' @export

compute_psd.eeg_epochs <- function(data,
                                   seg_length = NULL,
                                   noverlap = NULL,
                                   n_fft = 256,
                                   method = "Welch",
                                   keep_trials = TRUE,
                                   demean = TRUE,
                                   verbose = TRUE,
                                   ...) {

  if (demean) {
    data <- rm_baseline(data, verbose = verbose)
  }
  srate <- data$srate
  if (is.null(seg_length)) {
    seg_length <- n_fft
  }

  if (seg_length > n_fft) {
    stop("seg_length cannot be greater than n_fft")
  }
  data$signals <- split(data$signals,
                        data$timings$epoch)
  n_times <- nrow(data$signals[[1]])
  if (n_times < seg_length) {
    seg_length <- n_times
  }
  if (is.null(noverlap)) {
    noverlap <- seg_length %/% 2
  } else if (noverlap >= seg_length) {
    stop("noverlap should not be larger than seg_length.")
  }
  if (method == "Welch") {
    final_output <-
      lapply(data$signals,
             function(x) welch_fft(x,
                                   seg_length,
                                   noverlap = noverlap,
                                   n_fft = n_fft,
                                   srate = srate,
                                   n_sig = n_times)
             )
  }  else {
    stop("Welch is the only available method at this time.")
  }
  if (keep_trials) {
    final_output <- dplyr::bind_rows(final_output,
                                     .id = "epoch")
    final_output$epoch <- as.numeric(final_output$epoch)
    if (!is.null(epochs(data))) {
      final_output <- dplyr::left_join(final_output,
                                       epochs(data))
    }
  } else {
    final_output <- Reduce("+", final_output) / length(final_output)
  }
  final_output
}

#' @describeIn compute_psd Compute PSD for an \code{eeg_evoked} object
#' @export

compute_psd.eeg_evoked <- function(data,
                                   seg_length = NULL,
                                   noverlap = NULL,
                                   n_fft = 256,
                                   method = "Welch",
                                   demean = TRUE,
                                   verbose = TRUE,
                                   ...) {
  if (demean) {
    data <- rm_baseline(data, verbose = verbose)
  }
  srate <- data$srate

  if (is.null(seg_length)) {
    seg_length <- n_fft
  }

  if (seg_length > n_fft) {
    stop("seg_length cannot be greater than n_fft")
  }

  n_times <- nrow(data$signals)
  if (n_times < seg_length) {
    seg_length <- n_times
  }

  if (is.null(noverlap)) {
    noverlap <- seg_length %/% 8
  } else if (noverlap >= seg_length) {
    stop("noverlap should not be larger than seg_length.")
  }

    if (method == "Welch") {
    final_output <-
      welch_fft(data$signals,
                seg_length,
                noverlap = noverlap,
                n_fft = n_fft,
                srate = srate,
                n_sig = n_times)
  }  else {
    stop("Welch is the only available method at this time.")
  }

  final_output
}

#' Welch fft
#'
#' Internal function for calculating the PSD using Welch's method
#'
#' @param data Object to perform FFT on.
#' @param seg_length length of each segment of data.
#' @param n_fft length of FFT.
#' @param noverlap overlap between segments.
#' @param n_sig number of samples total.
#' @param srate Sampling rate of the data.
#' @importFrom stats fft
#' @keywords internal

welch_fft <- function(data,
                      seg_length,
                      n_fft,
                      noverlap,
                      n_sig,
                      srate) {

  # split data into segments
  if (seg_length < n_sig) {
    data_segs <- lapply(data,
                        split_vec,
                        seg_length,
                        noverlap)
    n_segs <- length(data_segs[[1]])
    #n_segs <- length(data_segs)
    # this splits the data into a list of ncol elements; each list element is
    # also a list containing n_segs elements - consider recoding this to combine
    # segments into
  } else {
    data_segs <- as.matrix(data)
    n_segs <- 1
  }

  # Hamming window.
  win <- .54 - (1 - .54) * cos(2 * pi * seq(0, 1, by = 1 / (seg_length - 1)))

  # Normalise the window
  U <- c(t(win) %*% win)

  #do windowing and zero padding if necessary, then FFT
  if (n_segs == 1) {
    data_segs <- sweep(data_segs,
                       1,
                       win, "*")

    if (n_fft > seg_length) {
       data_segs <- apply(data_segs, 2,
                          function(x) c(x,
                                        numeric(n_fft - seg_length)))
    }

    data_fft <- mvfft(data_segs)
    final_out <- apply(data_fft,
                       2,
                       function(x) abs(x * Conj(x)) / U)

    # Normalize by sampling rate
    if (is.null(srate)) {
      final_out <- final_out / (2 * pi)
      freqs <- seq(0, seg_length / 2) / (seg_length)
    } else {
      final_out <- final_out / srate
      freqs <- seq(0, n_fft / 2) / (n_fft) * srate
    }
  } else {
    data_segs <- lapply(data_segs,
                        function(x) lapply(x,
                                           function(y) y * win))
    if (n_fft > seg_length) {
      data_segs <- lapply(data_segs,
                          function(x) apply(data_segs,
                                            2,
                                            function(x) c(x,
                                                          numeric(n_fft - seg_length))))
    }
    data_fft <- lapply(data_segs,
                       function(x) lapply(x,
                                          fft))
    final_out <- lapply(data_fft,
                        function(x) sapply(x,
                                           function(y) abs(y * Conj(y)) / U))
    # Normalize by sampling rate or by signal length if no sampling rate
    if (is.null(srate)) {
      final_out <- rowMeans(as.data.frame(final_out)) / (2 * pi)
      freqs <- seq(0, seg_length / 2) / (seg_length)
    } else {
      final_out <- as.data.frame(lapply(final_out, rowMeans)) / srate
      freqs <- seq(0, n_fft / 2) / (n_fft) * srate
    }
  }

  #select first half of spectrum and double amps, output is power - uV^2 / Hz
  final_out <- final_out[1:(n_fft / 2 + 1), , drop = FALSE]
  #final_out[2:(n_fft / 2 + 1), ] <- (final_out[2:(n_fft / 2 + 1), ] * 2) ^ 2
  final_out <- data.frame(final_out,
                          frequency = freqs)
  final_out <- final_out[final_out$frequency > 0, ]
  final_out
}

#' Segment data
#'
#' Split data into segments for Welch PSD.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param vec Data vector to be split up into segments.
#' @param seg_length Length of segments to be FFT'd (in samples).
#' @param overlap Overlap between segments (in samples).
#' @keywords internal

split_vec <- function(vec, seg_length, overlap) {

  if (is.data.frame(vec)) {
    k <- floor((nrow(vec) - overlap) / (seg_length - overlap))
  } else {
    k <- floor((length(vec) - overlap) / (seg_length - overlap))
  }

  starts <- seq(1,
                k * (seg_length - overlap),
                by = seg_length - overlap)
  ends <- starts + seg_length - 1
  lapply(seq_along(starts),
         function(i) vec[starts[i]:ends[i]])
}




