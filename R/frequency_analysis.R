#' Compute power spectral density
#'
#' Returns PSD calculated using Welch's method for every channel in the data.
#' Currently concatenates across epochs. Returns uV^2 / Hz
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to be plotted. Accepts objects of class \code{eeg_data}
#' @param ... any further parameters passed to specific methods

compute_psd <- function(data, ...) {
  UseMethod("compute_psd", data)
}

#' @param n_fft Length of FFT to be calculated.
#' @param seg_length Length of rolling data segments.
#' @param noverlap Number of (sampling) points of overlap between segments.
#' @param srate Sampling rate
#' @param method Defaults to "Welch". No other method currently implemented.
#' @describeIn compute_psd Compute PSD for an \code{eeg_data} object

compute_psd.eeg_data <- function(data,
                                 seg_length = NULL,
                                 noverlap = 0,
                                 n_fft = 256,
                                 srate = NULL,
                                 method = "Welch") {

  srate <- data$srate

  if (is.null(seg_length)) {
    seg_length <- n_fft
  }

  if (seg_length > n_fft) {
    stop("seg_length cannot be greater than n_fft")
  }

  if (noverlap == 0) {
    noverlap <- seg_length %/% 8
  } else if (noverlap >= seg_length) {
    stop("noverlap should not be larger than seg_length.")
  }

  # split data into segments
  if (seg_length < nrow(data$signals)) {
    data_segs <- lapply(data$signals, split_vec, seg_length, noverlap)
    } else {
      data_segs <- data
    }

  # Hamming window.
  win <- .54 - (1 - .54) * cos(2 * pi * seq(0, 1, by = 1 / (n_fft - 1)))

  #do windowing and zero padding if necessary, then FFT
  if (n_fft > seg_length) {
    zero_pad <- rep(0, n_fft - seg_length)
    data_fft <- lapply(data_segs,
                       function(x) lapply(x,
                                          function(y) fft(c(y * win, zero_pad))))
  } else if (n_fft == seg_length){
    data_fft <- lapply(data_segs,
                       function(x) lapply(x,
                                          function(y) fft(y * win)))
  }

  # Normalise the window
  U <- c(t(win) %*% win)

  final_out <- lapply(data_fft,
                      function(x) sapply(x, function(y) abs(y * Conj(y)) / U))

  # Normalize by sampling rate
  if (is.null(srate)) {
    final_out <- rowMeans(as.data.frame(final_out)) / (2 * pi)
    freqs <- seq(0, seg_length / 2) / (seg_length)
    } else {
      final_out <- as.data.frame(lapply(final_out, rowMeans)) / srate
      freqs <- seq(0, n_fft / 2) / (n_fft) * srate
      }

  #select first half of spectrum and double amps, output is power - uV^2 / Hz
  final_out <- final_out[1:(n_fft / 2 + 1), ]
  final_out[2:(n_fft / 2 + 1), ] <- (final_out[2:(n_fft / 2 + 1), ] * 2) ^ 2
  data.frame(final_out, frequency = freqs)
}

#' Segment data.
#'
#' Split data into segments for Welch PSD.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param vec Data vector to be split up into segments.
#' @param seg_length Length of segments to be FFT'd (in samples).
#' @param overlap Overlap between segments (in samples).
#' @noRd

split_vec <- function(vec, seg_length, overlap) {
  k <- floor((length(vec) - overlap) / (seg_length - overlap))
  starts <- seq(1, k * (seg_length - overlap), by = seg_length - overlap)
  ends <- starts + seg_length - 1
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

#' Morlet wavelet
#'
#' Generate Morlet wavelet family
#'
#' @param frex Frequency range of interest
#' @param wavtime Times at which to place a wavelet
#' @param n_cycles Length of wavelet in cycles
#' @param n_freq number of frequencies to resolve
#' @noRd

morlet <- function(frex, wavtime, n_cycles = 7, n_freq = 30) {

  frex <- seq(frex[1], frex[2], length.out = n_freq)

  # widths of Gaussian
  if (length(n_cycles) == 1) {
    g_width <- n_cycles / (2 * pi * frex)
  } else {
    g_width <- seq(n_cycles[1], n_cycles[2], length.out = n_freq) / (2 * pi * frex)
  }

  t_by_f <- matrix(wavtime,
                   nrow = length(wavtime),
                   ncol = length(frex))

  c_sine <- 2 * 1i * pi * sweep(t_by_f, 2, frex, "*")
  gaussians <- sweep(t_by_f ^ 2, 2, 2 * g_width ^ 2, "/")
  m_family <- exp(c_sine - gaussians)
  m_family
}

#' N-point FFT
#'
#' Either zero-pads or truncates a signal to N and runs an FFT
#'
#' @param signal signal to be FFT'd
#' @param n Number of FFT points
#' @noRd

fft_n <- function(signal, n) {

  if (is.vector(signal)) {
    if (length(signal) < n) {
      signal_n <- c(signal, rep(0, n - length(signal)))
      } else {
        signal_n <- signal[1:n]
      }
    fft(signal_n)
  } else {
    if (nrow(signal) < n) {
      signal_n <- matrix(0, nrow = n, ncol = ncol(signal))
      signal_n[1:nrow(signal), ] <- signal
    } else {
      signal_n <- signal[1:n, ]
    }
    mvfft(as.matrix(signal_n))
  }

}

#' Convolve with morlets
#'
#' @param morlet_fam family of morlet wavelets
#' @param signal signal to be convolved
#' @param n points for FFT
#' @param wavtime time points
#' @param srate Sampling rate of the signal
#' @noRd

conv_mor <- function(morlet_fam, signal, n, wavtime, srate) {
  sigX <- fft_n(as.matrix(signal), n)

  tf_matrix <- array(dim = c(nrow(sigX), ncol(signal), ncol(morlet_fam)))

  for (i in 1:ncol(signal)) {
    tf_matrix[, i, ] <- mvfft(sigX[, i] * morlet_fam, inverse = TRUE) / srate
  }

  nHfkn <- floor(length(wavtime) / 2) + 1
  tf <- tf_matrix[nHfkn:(nrow(tf_matrix) - nHfkn + 1), , ]
  tf <- abs(tf) * 2
  tf
}

#' Time-frequency analysis
#'
#' Morlet wavelet time-frequency analysis.
#'
#' @param data EEG data to be TF transformed
#' @param ... Further parameters of the timefreq transformation
#' @export

tf_morlet <- function(data, ...) {
  UseMethod("tf_morlet", data)
}

#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param n_cycles Number of cycles at each frequency.
#' @param keep_trials Keep single trials or average over them before returning.
#' @importFrom abind abind
#' @describeIn tf_morlet Time-frequency decomposition of \code{eeg_epochs}
#'   object.
#' @export
tf_morlet.eeg_epochs <- function(data, foi, n_freq, n_cycles = 7, keep_trials = TRUE, ...) {

  if (length(foi) > 2) {
    stop("No more than two frequencies should be specified.")
  } else if (length(foi) == 2) {
    foi <- c(min(foi), max(foi))
  }

  wavtime <- unique(data$timings$time)

  morlet_family <- morlet(frex = foi,
                          n_freq = n_freq,
                          wavtime = wavtime,
                          n_cycles = n_cycles
                          )

  data$signals <- split(data$signals, data$timings$epoch)
  max_length <- nrow(data$signals[[1]])
  n_kern <- length(wavtime)
  n_conv <- max_length + n_kern - 1

  # zero-pad and run FFTs on morlets
  mf_zp <- fft_n(morlet_family, n_conv)

  # normalise wavelets
  mf_zp_maxes <- apply(abs(mf_zp), 2, which.max)
  mf_zp_maxes <- lapply(seq_along(mf_zp_maxes),
                  function(x) mf_zp[mf_zp_maxes[[x]], x])
  norm_mf <- matrix(unlist(lapply(seq_along(mf_zp_maxes),
                                  function(x) mf_zp[, x] / mf_zp_maxes[[x]])),
                    ncol = n_freq)

  data$signals <- lapply(data$signals,
                         function(x) conv_mor(norm_mf,
                                              x,
                                              n_conv,
                                              wavtime,
                                              data$srate))

  data$signals <- abind::abind(data$signals, along = 4)
  class(data) <- "eeg_tfr"
  data$freqs <- seq(foi[1], foi[2], length.out = n_freq)
  data

}
