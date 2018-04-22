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
    #data_segs <- split_vec(data$amplitude, seg_length, noverlap)
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

  final_out <- lapply(data_fft, function(x) sapply(x, function(y) abs(y * Conj(y)) / U))

  # Normalize by sampling rate
  if (is.null(srate)) {
    final_out <- rowMeans(as.data.frame(final_out)) / (2 * pi)
    freqs <- seq(0, seg_length / 2) / (seg_length)
    } else {
      final_out <- as.data.frame(lapply(final_out, rowMeans)) / srate
      freqs <- seq(0, n_fft / 2) / (n_fft) * srate
      }

  #select first half of spectrum and double amps, output is power - uV^2 / Hz
  final_out <- final_out[1:(n_fft / 2 + 1) , ]
  final_out[2:(n_fft / 2 + 1), ] <- (final_out[2:(n_fft / 2 + 1), ] * 2) ^2
  data.frame(final_out, frequency = freqs)
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

morlet <- function(frex, wavtime, n_cycles = 7, n_freq = 30) {

  frex <- seq(frex[1], frex[2], length.out = n_freq)

  # width of Gaussian
  if (length(n_cycles) == 1) {
    g_width <- n_cycles / (2 * pi * frex)
  } else {
    g_width <- seq(n_cycles[1], n_cycles[2], length.out = nFrex) / (2 * pi * frex)
  }

  t_by_f <- matrix(wavtime,
                   nrow = length(wavtime),
                   ncol = length(frex))
  c_sine <- 2* 1i * pi * t(t(t_by_f) * frex)
  gaussians <- t(t(-(t_by_f ^ 2)) / (2 * g_width ^ 2))
  m_family <- exp(c_sine + gaussians)
  # normalise wavelets
  m_max_index <- apply(abs(m_family), 2, which.max)
  m_max <- lapply(seq_along(m_max_index),
                  function(x) m_family[m_max_index[[x]], x])
  norm_mf <- matrix(unlist(lapply(seq_along(m_max),
                                  function(x) m_family[, x] / m_max[[x]])),
                    ncol = nFrex)
  norm_mf
}

#' N-point FFT
#'
#' Either zero-pads or truncates a signal to N and runs an FFT
#'
#' @param signal signal to be FFT'd
#' @param n Number of FFT points

fft_n <- function(signal, n) {
  if (length(signal) < n) {
    signal_n <- c(signal, rep(0, n - length(signal)))
  } else {
    signal_n <- signal[1:n]
  }
  fft(signal_n)
}

#' Convolve with morlets
#'
#' @param morlet_fam family of morlet wavelets
#' @param signal signal to be convolved
#' @param n points for FFT

conv_mor <- function(morlet_fam, signal, n) {
  sigX <- fft_n(signal, n)
  swe_test <- sweep(morlet_fam, 1, sigX, FUN = "*")
  sigdi <- mvfft(swe_test, inverse = TRUE)
  nHfkn <- floor(length(wavtime) / 2) + 1
  tf <- sigdi[nHfkn:(length(sigdi[, 1]) - nHfkn + 1), ]
  tf <- abs(tf) * 2
  tf
}


#' Time-frequency analysis
#'
#' Morlet wavelet time-frequency analysis.
#'
#' @param data EEG data to be TF transformed
#' @param ... Further parameters of the timefreq transformation

tf_morlet <- function(data, ...) {
  UseMethod("tf_morlet", data)
}

#' @param foi Frequencies of interest. Scalar or character vector of the lowest and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param n_cycles Number of cycles at each frequency.
#' @importFrom dplyr count
#'
#' @describeIn tf_morlet Time-frequency decomposition of \code{eeg_epochs} object.
tf_morlet.eeg_epochs <- function(data, foi, n_freq, n_cycles, ...) {

  if (length(foi) > 2) {
    stop("No more than two frequencies should be specified.")
  } else if (length(foi) == 2) {
    foi <- c(min(foi), max(foi))
  }

  morlet_family <- morlet(frex = foi,
                          n_freq = n_freq,
                          wavtime = unique(data$timings$time),
                          n_cycles = n_cycles
                          )
  max_length <- max(dplyr::count(eeg_epochs$timings, epoch)$n)

  # zero-pad before running ffts
  mf_zp <- apply(morlet_family, 2, function(x) c(x, rep(0, max_length - length(x))))
  morlet_fft <- mvfft(mf_zp)
  conv_mor(morlet_fft, signal, max_length)
}
