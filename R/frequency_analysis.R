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
#' points. If not specified, it defaults to 50% overlap between segments.
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
#' @describeIn compute_psd Compute PSD for an \code{eeg_data} object
#' @export

compute_psd.eeg_data <- function(data,
                                 seg_length = NULL,
                                 noverlap = NULL,
                                 n_fft = NULL,
                                 method = "Welch",
                                 ...) {

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
  if (is.null(noverlap)) {
    noverlap <- seg_length %/% 2
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

#' @describeIn compute_psd Compute PSD for an \code{eeg_evoked} object
#' @export

compute_psd.eeg_evoked <- function(data,
                                   seg_length = NULL,
                                   noverlap = NULL,
                                   n_fft = 256,
                                   method = "Welch",
                                   ...) {
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
  if (noverlap == 0) {
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
    n_segs <- length(data_segs)
    # this splits the data into a list of ncol elements; each list element is
    # also a list containing n_segs elements - consider recoding this to combine
    # segments into
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
    data_segs <- sweep(data_segs, 1, win, "*")

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
  final_out[2:(n_fft / 2 + 1), ] <- (final_out[2:(n_fft / 2 + 1), ] * 2) ^ 2
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


#' Compute Time-Frequency representation of EEG data
#'
#' This function creates a time frequency represention of EEG time series data.
#' Currently, the only available method is a Morlet wavelet transformation
#' performed using convolution in the frequency domain.
#'
#' @param data An object of class \code{eeg_epochs}.
#' @param ... Further TFR parameters
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @examples
#' compute_tfr(demo_epochs, method = "morlet", foi = c(4, 30), n_freq = 10, n_cycles = 3)
#' @export

compute_tfr <- function(data, ...) {
  UseMethod("compute_tfr", data)
}

#' @describeIn compute_tfr Default method for compute_tfr
#' @export
compute_tfr.default <- function(data, ...) {
  warning("compute_tfr requires data in eeg_epochs format.")
}

#' @param method Time-frequency analysis method. Defaults to "morlet".
#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved. Scalar.
#' @param n_cycles Scalar. Number of cycles at each frequency. Currently only
#'   supports a single number of cycles at all frequencies.
#' @param keep_trials Keep single trials or average over them before returning.
#'   Defaults to FALSE.
#' @param output Sets whether output is power, phase, or fourier coefficients.
#' @param downsample Downsampling factor. Integer. Selects every n samples after
#'   performing time-frequency analysis.
#' @describeIn compute_tfr Default method for compute_tfr
#' @export

compute_tfr.eeg_epochs <- function(data,
                                   method = "morlet",
                                   foi,
                                   n_freq,
                                   n_cycles = 7,
                                   keep_trials = FALSE,
                                   output = "power",
                                   downsample = 1,
                                   ...) {


  tfr_obj <- switch(method,
                   "morlet" = tf_morlet(data,
                                        foi,
                                        n_freq,
                                        n_cycles,
                                        keep_trials,
                                        output = output,
                                        downsample),
                   warning("Unknown method supplied. Currently supported method is 'morlet'"))
  tfr_obj
}

#' @describeIn compute_tfr Method for \code{eeg_evoked} objects.
compute_tfr.eeg_evoked <- function(data,
                                   method = "morlet",
                                   foi,
                                   n_freq,
                                   n_cycles = 7,
                                   keep_trials = FALSE,
                                   output = "power",
                                   downsample = 1,
                                   ...) {
  switch(method,
         "morlet" = tf_morlet(data,
                              foi,
                              n_freq,
                              n_cycles,
                              keep_trials,
                              output,
                              downsample),
         warning("Unknown method supplied. Currently supported method is 'morlet'"))
}


#' Perform Morlet time-frequency analysis
#'
#' Internal function for performing Morlet wavelet transforms using convolution
#' in frequency domain
#'
#' @param data Data in \code{eeg_epochs} format.
#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param n_cycles Number of cycles at each frequency.
#' @param keep_trials Keep single trials or average over them before returning.
#' @param output Sets whether output is power, phase, or fourier coefficients.
#' @param downsample Downsampling factor (integer)
#' @importFrom abind abind
#' @keywords internal

tf_morlet <- function(data,
                      foi,
                      n_freq,
                      n_cycles,
                      keep_trials,
                      output,
                      downsample) {

  if (length(foi) > 2) {
    stop("No more than two frequencies should be specified.")
  } else if (length(foi) == 2) {
    foi <- c(min(foi),
             max(foi))
  } else {
    foi <- c(foi, foi)
    n_freq <- 1
  }

  frex <- seq(foi[1],
              foi[2],
              length.out = n_freq)

  # if a min and max n_cycles is specified, expand out to cycles per n_freq
  if (length(n_cycles) == 2) {
    n_cycles <- seq(n_cycles[1],
                    n_cycles[2],
                    length.out = n_freq)
  } else if (length(n_cycles) > 2) {
    stop("n_cycles should be a vector of length 1 or length 2.")
  }

  #de-mean each epoch
  #data <- rm_baseline(data)
  elecs <- names(data$signals)
  fft_points <- length(unique(data$timings$time))
  sigtime <- unique(data$timings$time)

  #Create a family of morlet wavelets (unscaled)
  morlet_family <- morlet(frex = frex,
                          srate = data$srate,
                          n_freq = n_freq,
                          n_cycles = n_cycles
                          )

  # This is a total hack to make the rest of the code behave with eeg_evoked data
  if (is.eeg_evoked(data)) {
    data$timings$epoch <- 1
  }

  data$signals <- split(data$signals,
                        data$timings$epoch)
  max_length <- nrow(data$signals[[1]])
  n_kern <- nrow(morlet_family)
  n_conv <- max_length + n_kern - 1

  # zero-pad and run FFTs on morlets
  mf_zp <- fft_n(morlet_family,
                 n_conv)


  # Normalise wavelets for FFT (as suggested by Mike X. Cohen):
  norm_mf <- wavelet_norm(mf_zp,
                          n_freq)

  # Run the FFT convolutions on each individual trial

  # generate a vector for selection of specific timepoints
  time_sel <- seq(1, max_length, by = downsample)

  trial_conv <- function(trial_dat) {
    trial_dat <-
      convert_tfr(
        conv_mor(norm_mf,
                 trial_dat,
                 n_conv,
                 sigtime,
                 data$srate),
        output)
    trial_dat[time_sel, , , drop = FALSE]
  }

  data$signals <- lapply(data$signals,
                         trial_conv)

  sigtime <- sigtime[time_sel]
  data$timings <- data$timings[data$timings$time %in% sigtime, ]
  n_epochs <- length(data$signals)
  sig_dims <- dim(data$signals[[1]])

  # Bind single trial matrices together into a single matrix
  # do.call method is slower than abind, but uses less memory
  data$signals <- do.call(rbind, data$signals)

  dim(data$signals) <- c(n_epochs, sig_dims)

  # Create a list of metadata about the TFR
  data$freq_info <- list(freqs = frex,
                         morlet_resolution = morlet_res(frex,
                                                        n_cycles),
                         method = "morlet",
                         output = output,
                         baseline = "none")

  # Remove edges of the TFR'd data, where edge effects would be expected.
  edge_mat <- remove_edges(sigtime,
                           data$freq_info$morlet_resolution$sigma_t)

  if (keep_trials) {
    data$signals <- aperm(data$signals,
                          c(2, 3, 4, 1))
    dimnames(data$signals) <- list(sigtime,
                                   elecs,
                                   data$freq_info$freqs,
                                   unique(data$timings$epoch))
    data$signals <- sweep(data$signals,
                          c(1, 3),
                          edge_mat,
                          "*")
    data <- eeg_tfr(data$signals,
                    srate = data$srate,
                    events = data$events,
                    chan_info = data$chan_info,
                    reference = data$reference,
                    timings = data$timings,
                    freq_info = data$freq_info,
                    dimensions = c("time",
                                   "electrode",
                                   "frequency",
                                   "epoch"))
    return(data)
  }


  data <- average_tf(data)

  dimnames(data$signals) <- list(sigtime,
                                 elecs,
                                 data$freq_info$freqs)
  data$signals <- sweep(data$signals,
                        c(1, 3),
                        edge_mat,
                        "*")
  data <- eeg_tfr(data$signals,
                  srate = data$srate,
                  events = NULL,
                  chan_info = data$chan_info,
                  reference = data$reference,
                  timings = data$timings[1:nrow(data$signals), "time"],
                  freq_info = data$freq_info,
                  dimensions = c("time",
                                 "electrode",
                                 "frequency"))
  data
}

#' Morlet wavelet
#'
#' Generate Morlet wavelet family
#'
#' @param frex Frequency range of interest
#' @param n_cycles Length of wavelet in cycles
#' @param n_freq number of frequencies to resolve
#' @keywords internal

morlet <- function(frex,
                   srate,
                   n_cycles,
                   n_freq) {

  # calculate frequency and temporal std devs
  sigma_t <- morlet_res(frex, n_cycles)$sigma_t
  max_sd_t <- max(sigma_t)
  tstep <- 1 / srate

  # round the max SD to the next biggest number divisible by tstep
  round_sd <- max_sd_t + (tstep - (max_sd_t %% tstep))

  # calculate length of kernel as 6 * maximum temporal SD
  # TO DO - change this to do it for every frequency separately
  wavtime <- seq(-round_sd * 3,
                 round_sd * 3,
                 by = tstep)

  t_by_f <- matrix(wavtime,
                   nrow = length(wavtime),
                   ncol = length(frex))

  # Create sine waves at each frequency
  c_sine <- 2 * 1i * pi * sweep(t_by_f,
                                2,
                                frex,
                                "*")
  gaussians <- sweep(t_by_f ^ 2,
                     2,
                     2 * sigma_t ^ 2,
                     "/")
  m_family <- exp(c_sine - gaussians)
 # A <- round_sd^(-1/2) * pi^(-1/4) # Equivalent to Fieldtrip scaling (1/sqrt(sd * sqrt(pi)))
  m_family
}

#' Morlet wavelet resolutions
#'
#' Calculate frequency and temporal standard deviations for the Morlet wavelets
#' @param frex Frequencies of interest
#' @param n_cycles Number of cycles for each frequency
#' @keywords internal
morlet_res <- function(frex, n_cycles) {
  sigma_f <- frex / n_cycles
  sigma_t <- 1 / (2 * pi * sigma_f)
  data.frame(sigma_f, sigma_t)
}

#' N-point FFT
#'
#' Either zero-pads by adding trailing zeroes, or truncates a signal to N and
#' runs an FFT.
#'
#' @param signal signal to be FFT'd
#' @param n Number of FFT points
#' @keywords internal

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
      signal_n <- matrix(0,
                         nrow = n,
                         ncol = ncol(signal))
      signal_n[1:nrow(signal), ] <- signal
    } else {
      signal_n <- signal[1:n, ]
    }
    mvfft(as.matrix(signal_n))
  }

}

#' Convolve with morlets
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param morlet_fam family of morlet wavelets
#' @param signal signal to be convolved
#' @param n points for FFT
#' @param wavtime time points
#' @param srate Sampling rate of the signal
#' @importFrom stats mvfft
#' @keywords internal

conv_mor <- function(morlet_fam,
                     signal,
                     n,
                     wavtime,
                     srate) {

  n_chans <- ncol(signal)
  n_times <- nrow(signal)
  signal <- fft_n(as.matrix(signal),
                n)
  tf_matrix <- array(1i, dim = c(nrow(signal),
                                 n_chans,
                                 ncol(morlet_fam)))
  for (i in 1:ncol(signal)) {
    tf_matrix[, i, ] <- mvfft(signal[, i] * morlet_fam,
                              inverse = TRUE) / srate
  }

  nkern <- n - n_times + 1
  nHfkn <- floor(nkern / 2) + 1
  tf_matrix <- tf_matrix[nHfkn:(nrow(tf_matrix) - nHfkn + 1), , , drop = FALSE]
  tf_matrix
}



#' Remove convolution edges
#'
#' Create a matrix indicating which timepoints likely suffer from edge effects.
#' Returns a time by frequency matrix with NA
#' @param sigtime timepoints in the signal
#' @param sigma_t standard deviations of the morlet wavelets
#' @keywords internal

remove_edges <- function(sigtime, sigma_t) {
  #Full width of wavelet kernel is 2 * 3 * sigma_t
  sigma_t <- sigma_t * 3
  #create a matrix of timepoints for every frequency
  times_mat <- matrix(sigtime,
                      nrow = length(sigtime),
                      ncol = length(sigma_t))
  #calculate time of left edge and replace anything earlier with NA
  left_edge <- sigtime[[1]] + sigma_t
  times_mat <- t(apply(times_mat,
                       1,
                       function(x) ifelse(x < left_edge,
                                          NA,
                                          x)))
  #calculate time of right edge and replace anything later with NA
  right_edge <- sigtime[[length(sigtime)]] - sigma_t
  times_mat <- t(apply(times_mat,
                       1,
                       function(x) ifelse(x > right_edge,
                                          NA,
                                          x)))
  #Finally, replace anything that isn't NA with 1
  times_mat <- ifelse(is.na(times_mat), NA, 1)
  times_mat
}


#' Calculate circular mean
#'
#' Calculates the circular mean from vector of phase angles. Input should be in
#' radians.
#'
#' @param data vector of phase angles in radians
#' @keywords internal
circ_mean <- function(data) {
  atan(mean(sin(data)) / mean(cos(data)))
}

#' Convert Fourier output to power or phase as requested.
#'
#' @param data Fourier coefficients from from eeg_tfr
#' @param output What output is desired - "power", "phase" or "fourier"
#' @keywords internal
convert_tfr <- function(data, output) {
  if (is.complex(data[[1]])) {
    switch(output,
         "power" = abs(data) ^ 2,
         "phase" = atan2(Im(data),
                         Re(data)),
         "fourier" = data)
  } else {
    warning("Data is not complex, returning original data.")
  }
}

#' Normalise Morlet wavelets
#'
#' Normalise the FFT'd Morlet wavelet family as recommended by Mike X Cohen,
#' dividing each wavelet by its absolute maximum. This should result in each
#' frequency being passed at unit amplitude, and the resulting convolution with
#' the signal should return units approximately on the original scale (i.e. uV^2
#' / Hz)
#'
#' @param mf_zp A zero-padded, FFT'd morlet wavelet family
#' @param n_freq Number of frequencies
#' @keywords internal
wavelet_norm <- function(mf_zp, n_freq) {
  mf_zp_maxes <- apply(abs(mf_zp),
                       2,
                       which.max)
  mf_zp_maxes <- lapply(seq_along(mf_zp_maxes),
                        function(x) mf_zp[mf_zp_maxes[[x]], x])
  norm_mf <- lapply(seq_along(mf_zp_maxes),
                    function(x) mf_zp[, x] / mf_zp_maxes[[x]])
  norm_mf <- matrix(unlist(norm_mf),
                    ncol = n_freq)
}

#' Internal function for averaging over epochs
#' @param data data to average over
#' @keywords internal
average_tf <- function(data) {

  if (data$freq_info$output == "phase") {
    data$signals <- apply(data$signals,
                          c(1, 2, 3),
                          circ_mean)
  } else {
    avg_tf <- array(0, dim = dim(data$signals)[2:4])
    for (iz in 1:dim(data$signals)[3]) {
      for (ij in 1:dim(data$signals)[4]) {
        avg_tf[, iz, ij] <- colMeans(data$signals[ , , iz, ij, drop = FALSE])
      }
    }
    data$signals <- avg_tf
  }
  data
}
