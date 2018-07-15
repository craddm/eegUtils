#' Compute power spectral density
#'
#' Returns PSD calculated using Welch's method for every channel in the data.
#' Currently concatenates across epochs. Returns uV^2 / Hz
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to be plotted. Accepts objects of class \code{eeg_data}
#' @param ... any further parameters passed to specific methods
#' @export

compute_psd <- function(data, ...) {
  UseMethod("compute_psd", data)
}

#' @param n_fft Length of FFT to be calculated.
#' @param seg_length Length of rolling data segments.
#' @param noverlap Number of (sampling) points of overlap between segments.
#' @param srate Sampling rate
#' @param method Defaults to "Welch". No other method currently implemented.
#' @describeIn compute_psd Compute PSD for an \code{eeg_data} object
#' @export

compute_psd.eeg_data <- function(data,
                                 seg_length = NULL,
                                 noverlap = 0,
                                 n_fft = 256,
                                 srate = NULL,
                                 method = "Welch",
                                 ...) {

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

  if (method == "Welch") {
    final_output <- #lapply(data$signals, function(x)
      welch_fft(data$signals,
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

#' @param keep_trials Include FFT for every trial in output
#' @describeIn compute_psd Compute PSD for an \code{eeg_epochs} object
#' @export

compute_psd.eeg_epochs <- function(data,
                                   seg_length = NULL,
                                   noverlap = 0,
                                   n_fft = 256,
                                   srate = NULL,
                                   method = "Welch",
                                   keep_trials = TRUE,
                                   ...) {
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
  if (noverlap == 0) {
    noverlap <- seg_length %/% 8
  } else if (noverlap >= seg_length) {
    stop("noverlap should not be larger than seg_length.")
  }
  if (method == "Welch") {
    final_output <- lapply(data$signals, function(x)
      welch_fft(x,
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
    final_output <- dplyr::bind_rows(final_output, .id = "epoch")
  } else {
    final_output <- Reduce("+", final_output) / length(final_output)
    final_output
  }
}

#' Welch fft
#'
#' @param data Object to perform FFT on
#' @param seg_length length of each segment of data
#' @param n_fft length of FFT
#' @param noverlap overlap between segments
#' @param n_sig number of samples total
#' @param srate Sampling rate of the data
#' @importFrom stats fft
#' @noRd

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
    n_segs <- length(data_segs)
    # this splits the data into a list of ncol elements; each list element is
    # also a list containing n_segs elements - consider recoding this to combine segments into
  } else {
    data_segs <- data
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
                       win,
                       "*")
    if (n_fft > seg_length) {
      data_segs <- apply(data_segs,
                         2,
                         function(x) c(x,
                                       numeric(n_fft - seg_length)))
    }
    data_fft <- mvfft(data_segs)
    final_out <- apply(data_fft,
                       2,
                       function(x) abs(x * Conj(x)) / U)
    # Normalize by sampling rate
    if (is.null(srate)) {
      final_out <- sweep(final_out,
                         1,
                         (2 * pi),
                         "/")
      freqs <- seq(0, seg_length / 2) / (seg_length)
    } else {
      #final_out <- as.data.frame(final_out) / srate
      final_out <- sweep(final_out,
                         1,
                         srate,
                         "/")
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
    # Normalize by sampling rate
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
  final_out[2:(n_fft / 2 + 1), ] <- (final_out[2:(n_fft / 2 + 1), ] * 2) ^ 2
  final_out <- data.frame(final_out,
                          frequency = freqs)
  final_out <- final_out[final_out$frequency > 0, ]
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
  if (is.data.frame(vec)) {
    k <- floor((nrow(vec) - overlap) / (seg_length - overlap))
    # starts <- seq(1, k * (seg_length - overlap), by = seg_length - overlap)
    # ends <- starts + seg_length - 1
    # lapply(seq_along(starts), function(i) vec[starts[i]:ends[i]])
  } else {
    k <- floor((length(vec) - overlap) / (seg_length - overlap))
  }
    starts <- seq(1,
                  k * (seg_length - overlap),
                  by = seg_length - overlap)
    ends <- starts + seg_length - 1
    lapply(seq_along(starts), function(i) vec[starts[i]:ends[i]])
}


#' Compute Time-Frequency representation of EEG data
#'
#' This function creates a time frequency represention of EEG time series data.
#' Currently, the only available method is a Morlet wavelet transformation
#' performed using convolution in the frequency domain.
#'
#' @param data An object of class \code{eeg_epochs}.
#' @param ... Furthere TFR parameters
#' @author Matt Craddock \email{matt@@mattcraddock.com}

compute_tfr <- function(data, ...) {
  UseMethod("compute_tfr", data)
}

#' @describeIn compute_tfr Default method for compute_tfr

compute_tfr.default <- function(data, ...) {
  warning("compute_tfr requires data in eeg_epochs format.")
}

#' @param method Time-frequency analysis method. Defaults to "morlet".
#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param n_cycles Number of cycles at each frequency.
#' @param keep_trials Keep single trials or average over them before returning.
#' @describeIn compute_tfr Default method for compute_tfr

compute_tfr.eeg_epochs <- function(data,
                                   method = "morlet",
                                   foi,
                                   n_freq,
                                   n_cycles = 7,
                                   keep_trials = TRUE, ...) {

  switch(method,
         "morlet" = tf_morlet(data,
                              foi,
                              n_freq,
                              n_cycles,
                              keep_trials),
         warning("Unknown method supplied. Currently supported method is 'morlet'"))

}

#' @param data Data in \code{eeg_epochs} format.
#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param n_cycles Number of cycles at each frequency.
#' @param keep_trials Keep single trials or average over them before returning.
#' @importFrom abind abind
#' @noRd
tf_morlet <- function(data,
                      foi,
                      n_freq = 10,
                      n_cycles = 7,
                      keep_trials = TRUE) {

  if (length(foi) > 2) {
    stop("No more than two frequencies should be specified.")
  } else if (length(foi) == 2) {
    foi <- c(min(foi),
             max(foi))
  }

  data <- rm_baseline(data)
  elecs <- names(data$signals)
  frex <- seq(foi[1], foi[2], length.out = n_freq)
  mod_frex <- frex %% 1
  freq_res <- min(mod_frex[mod_frex > 0])
  fft_points <- ceiling(data$srate / freq_res)
  wavtime <- unique(data$timings$time)
  sigtime <- unique(data$timings$time)
  morlet_family <- morlet(frex = foi,
                          srate = data$srate,
                          n_freq = n_freq,
                          wavtime = wavtime,
                          n_cycles = n_cycles
                          )

  data$signals <- split(data$signals,
                        data$timings$epoch)
  max_length <- nrow(data$signals[[1]])
  n_kern <- fft_points
  n_conv <- max_length + n_kern - 1
  #if ((n_conv / 2) %% 1 > 0) {
   # ceiling((max_length + n_kern - 1) / 2) * 2

  # zero-pad and run FFTs on morlets
  mf_zp <- fft_n(morlet_family,
                 n_conv)

  # Normalise wavelets:
  # 1) get the index for the absolute maximum
  # 2) divide each wavelet by its absolute maximum
  mf_zp_maxes <- apply(abs(mf_zp),
                       2,
                       which.max)
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

  data$signals <- abind::abind(data$signals,
                               along = length(dim(data$signals[[1]])) + 1)
  class(data) <- "eeg_tfr"
  data$freqs <- seq(foi[1],
                    foi[2],
                    length.out = n_freq)
  if (keep_trials) {
    dimnames(data$signals) <- list(sigtime,
                                   elecs,
                                   data$freqs,
                                   unique(data$timings$epoch))
    data
  } else {
      data$signals <- apply(data$signals,
                            c(1, 2, 3),
                            mean)
      dimnames(data$signals) <- list(sigtime,
                                     elecs,
                                     data$freqs)
      data
    }
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

morlet <- function(frex,
                   wavtime,
                   srate,
                   n_cycles,
                   n_freq) {

  frex <- seq(frex[1],
              frex[2],
              length.out = n_freq)
  frac_freq <- (frex %% 1)
  required_res <- min(frac_freq[frac_freq > 0])
  n_points <- srate / required_res

  orig_wavtime <- wavtime
  #n_points <- length(orig_wavtime)
  #n_points <- 512
  wavtime <- seq(0, 1 / srate * n_points, by = 1 / srate)
  wavtime <- wavtime - stats::median(wavtime)
  # widths of Gaussian
  if (length(n_cycles) == 1) {
    g_width <- n_cycles / (2 * pi * frex)
  } else {
    g_width <- seq(n_cycles[1],
                   n_cycles[2],
                   length.out = n_freq) / (2 * pi * frex)
  }

  t_by_f <- matrix(wavtime,
                   nrow = length(wavtime),
                   ncol = length(frex))
  c_sine <- 2 * 1i * pi * sweep(t_by_f,
                                2,
                                frex,
                                "*")
  gaussians <- sweep(t_by_f ^ 2,
                     2,
                     2 * g_width ^ 2,
                     "/")
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
#' @importFrom stats mvfft
#' @noRd

conv_mor <- function(morlet_fam,
                     signal,
                     n,
                     wavtime,
                     srate) {

  sigX <- fft_n(as.matrix(signal),
                n)
  tf_matrix <- array(dim = c(nrow(sigX),
                             ncol(signal),
                             ncol(morlet_fam)))
  for (i in 1:ncol(signal)) {
    tf_matrix[, i, ] <- mvfft(sigX[, i] * morlet_fam, inverse = TRUE) / srate
  }

  nkern <- n - nrow(signal) + 1
  nHfkn <- floor(nkern / 2) + 1
  tf <- tf_matrix[nHfkn:(nrow(tf_matrix) - nHfkn + 1), , ]
  tf <- abs(tf) * 2
  tf
}
