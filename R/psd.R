#' Compute power spectral density using Welch's method
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param df Data to be plotted. Accepts objects of class \code{eeg_data} or
#'   simple vectors. If a vector is supplied, \code{srate} must also be
#'   provided.
#' @param n_fft Length of FFT to be calculated.
#' @param seg_length Length of rolling data segments.
#' @param noverlap Number of (sampling) points of overlap between segments.
#' @param srate Sampling rate
#' @import tidyr
#' @import dplyr
#' @importFrom purrr map
#'

compute_psd_welch <- function(data,
                              seg_length = NULL,
                              noverlap = 0,
                              n_fft = 256,
                              srate = NULL) {

  if (is.eeg_data(data)) {
    srate <- data$srate
    data <- as.data.frame(data, long = TRUE)
  }

  if (is.null(seg_length)) {
    seg_length <- n_fft
  }

  if (seg_length < length(data)) {
    data_segs <- split_vec(data, seg_length, noverlap)
    } else {
      data_segs <- data
      }

  # Hamming window.
  win <- .54 - (1 - .54) * cos(2 * pi * seq(0, 1, by = 1 / (n_fft - 1)))

  if (seg_length > n_fft) {
    error('Seg_length must be >= n_fft')
    } else if (n_fft > seg_length) {
      zero_pad <- rep(0, n_fft)
      data_segs <- purrr::map(data_segs, ~ zero_pad + (. * win))
      data_fft <- purrr::map(data_segs, ~ fft(.))
    } else {
      data_fft <- purrr::map(data_segs, ~ fft(. * win))
    }

  # Normalise the window
  U <- c(t(win) %*% win)

  final_out <- purrr::map(data_fft, ~abs(. * Conj(.)) / U)

  if (is.null(srate)) {
    final_out <- rowMeans(as.data.frame(final_out)) / (2*pi)
    freqs <- seq(0, seg_length/2) / (seg_length)
    } else {
      final_out <- rowMeans(as.data.frame(final_out)) / srate
      freqs <- seq(0, n_fft/2) / (n_fft) * srate
      }

  final_out <- final_out[1:(n_fft/2 +1)]
  final_out[2:(n_fft/2 + 1)] <- final_out[2:(n_fft/2 +1)] * 2
  data.frame(power = final_out, frequency = freqs)
}


#' Segment data.
#'
#' Split data into segments for Welch PSD.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param vec Data vector to be split up into segments.
#' @param seg_length Length of segments to be FFT'd (in samples).
#' @param overlap Overlap between segments (in samples).

split_vec <- function(vec, seg_length, overlap) {
  k <- floor((length(vec) - overlap) / (seg.length - overlap))
  starts <- seq(1, k * (seg.length-overlap), by = seg.length-overlap)
  ends <- starts + seg.length - 1
  #ends[ends > length(vec)] <- length(vec)
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}
