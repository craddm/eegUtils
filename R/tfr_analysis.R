#' Compute Time-Frequency representation of EEG data
#'
#' This function creates a time frequency representation of EEG time series
#' data. Currently, it is possible to use either Morlet wavelets or Hanning
#' tapers during the decomposition, which uses convolution in the frequency
#' domain.
#'
#' @param data An object of class `eeg_epochs`.
#' @param ... Further TFR parameters
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @return An object of class `eeg_tfr`
#' @examples
#' out <- compute_tfr(demo_epochs, method = "morlet", foi = c(4, 30), n_freq = 10, n_cycles = 3)
#' out
#' out$freq_info$morlet_resolution
#' out <- compute_tfr(demo_epochs, method = "morlet", foi = c(4, 30), n_freq = 10, n_cycles = c(3, 10))
#' out$freq_info$morlet_resolution
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
#' @param n_freq Number of frequencies to be resolved. Must be an integer number
#'   of frequencies.
#' @param n_cycles Number of cycles at each frequency. If a single integer, use
#'   a constant number of cycles at each frequency. If a character vector of
#'   length 2, the number of cycles will scale with frequency from the minimum
#'   to the maximum.
#' @param spacing Use "linear" or "log" spacing for the frequency vector and number of cycles.
#' @param keep_trials Keep single trials or average over them before returning.
#'   Defaults to FALSE.
#' @param output Sets whether output is power, phase, or fourier coefficients.
#' @param downsample Downsampling factor. Integer. Selects every n samples after
#'   performing time-frequency analysis on the full sampling rate data.
#' @param verbose Print informative messages in console.
#' @describeIn compute_tfr Default method for `compute_tfr`
#' @export

compute_tfr.eeg_epochs <- function(data,
                                   method = "morlet",
                                   foi,
                                   n_freq,
                                   spacing = "linear",
                                   n_cycles = 7,
                                   keep_trials = FALSE,
                                   output = "power",
                                   downsample = 1,
                                   verbose = TRUE,
                                   ...) {
  if (identical(output, "fourier") && keep_trials == FALSE) {
    if (verbose) {
      message("For fourier output, all trials are kept.")
    }
    keep_trials <- TRUE
  }

  tfr_obj <- switch(
    method,
    "morlet" = tf_morlet(
      data = data,
      foi = foi,
      n_freq = n_freq,
      spacing = spacing,
      n_cycles = n_cycles,
      keep_trials = keep_trials,
      output = output,
      downsample = downsample,
      verbose = verbose
      ),
    "hanning" = tf_hanning(
      data = data,
      foi = foi,
      n_freq = n_freq,
      spacing = spacing,
      n_cycles = n_cycles,
      keep_trials = keep_trials,
      output = output,
      downsample = downsample,
      verbose = verbose
    ),
    warning(
      "Unknown method supplied. Currently supported methods are 'morlet', and 'hanning'"
    )
  )
  tfr_obj
}

#' @describeIn compute_tfr Method for `eeg_evoked` objects.
#' @export
compute_tfr.eeg_evoked <- function(data,
                                   method = "morlet",
                                   foi,
                                   n_freq,
                                   spacing = "linear",
                                   n_cycles = 7,
                                   keep_trials = FALSE,
                                   output = "power",
                                   downsample = 1,
                                   verbose = TRUE,
                                   ...) {
  if (identical(output, "fourier") && keep_trials == FALSE) {
    if (verbose) {
      message("For fourier output, all trials are kept.")
    }
    keep_trials <- TRUE
  }
  switch(
    method,
    "morlet" = tf_morlet(
      data,
      foi = foi,
      n_freq = n_freq,
      spacing = spacing,
      n_cycles = n_cycles,
      keep_trials = keep_trials,
      output = output,
      downsample = downsample,
      verbose = verbose
    ),
    "hanning" = tf_hanning(
      data = data,
      foi = foi,
      n_freq = n_freq,
      spacing = spacing,
      n_cycles = n_cycles,
      keep_trials = keep_trials,
      output = output,
      downsample = downsample,
      verbose = verbose
    ),
    warning(
      "Unknown method supplied. Currently supported methods are 'morlet' and 'hanning'"
    )
  )
}

#' @export
compute_tfr.eeg_group <- function(data,
                                  method = "morlet",
                                  foi,
                                  n_freq,
                                  n_cycles = 7,
                                  keep_trials = FALSE,
                                  output = "power",
                                  downsample = 1,
                                  verbose = TRUE,
                                  ...) {
  stop("Cannot compute tfr for grouped data.")
}

#' Perform Morlet time-frequency analysis
#'
#' Internal function for performing Morlet wavelet transforms using convolution
#' in frequency domain
#'
#' @param data Data in `eeg_epochs` format.
#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param spacing Use linear or log spacing for frequencies.
#' @param n_cycles Number of cycles at each frequency.
#' @param keep_trials Keep single trials or average over them before returning.
#' @param output Sets whether output is power, phase, or fourier coefficients.
#' @param downsample Downsampling factor (integer).
#' @param demean Remove mean before transforming.
#' @param verbose Print informative messages in console.
#' @importFrom abind abind
#' @importFrom stats nextn
#' @keywords internal

tf_morlet <- function(data,
                      foi,
                      n_freq,
                      spacing,
                      n_cycles,
                      keep_trials,
                      output,
                      downsample,
                      demean = TRUE,
                      verbose) {

  if (verbose) {
    message("Computing TFR using Morlet wavelet convolution")
  }

  frex <- parse_frex(foi,
                     n_freq,
                     verbose,
                     spacing = spacing)

  if (downsample > 1 && verbose) {
    message("Downsampling after transformation by a factor of ", downsample)
  }

  n_freq <- length(unique(frex))
  # if a min and max n_cycles is specified, expand out to cycles per n_freq
  if (length(n_cycles) == 2) {
    if (identical(spacing, "linear")) {
      n_cycles <- seq(n_cycles[1],
                      n_cycles[2],
                      length.out = n_freq)
    } else if (identical(spacing, "log")) {
      n_cycles <- exp(
        seq(log(n_cycles[1]),
            log(n_cycles[2]),
            length.out = n_freq)
        )
    }
      } else if (length(n_cycles) > 2) {
    stop("n_cycles should be a vector of length 1 or length 2.")
  }

  #de-mean each epoch
  if (demean) {
    data <- rm_baseline(data,
                        verbose = verbose)
  }

  elecs <- names(data$signals)
  #fft_points <- length(unique(data$timings$time))
  sigtime <- unique(data$timings$time)

  #Create a family of morlet wavelets (unscaled)
  morlet_family <- morlet(
    frex = frex,
    srate = data$srate,
    n_freq = n_freq,
    n_cycles = n_cycles
  )

  # Create a list of metadata about the TFR
  data$freq_info <- list(
    freqs = frex,
    morlet_resolution = morlet_res(frex,
                                   n_cycles),
    method = "morlet",
    output = output,
    baseline = "none"
  )

  n_kern <- nrow(morlet_family)
  max_length <- length(unique(data$timings$time))
  n_conv <- max_length + n_kern - 1
  n_conv <- stats::nextn(n_conv, 2)

  # zero-pad and run FFTs on morlets
  mf_zp <- fft_n(morlet_family,
                 n_conv)

  # Normalise wavelets for FFT (as suggested by Mike X. Cohen):
  norm_mf <- wavelet_norm(mf_zp,
                          n_freq)

  # Run the FFT convolutions on each individual trial

  # generate a vector for selection of specific timepoints
  time_sel <- seq(1,
                  length(unique(data$timings$time)),
                  by = downsample)

  data$signals <- run_tf(data,
                         norm_mf,
                         n_conv,
                         n_kern,
                         keep_trials,
                         time_sel,
                         output)

  sigtime <- sigtime[time_sel]
  data$timings <- data$timings[data$timings$time %in% sigtime,]

  if (keep_trials == TRUE) {
    dimnames(data$signals) <- list(
      epoch = unique(data$timings$epoch),
      time = sigtime,
      electrode = elecs,
      frequency = data$freq_info$freqs
    )
  } else {
    data$signals <- array(data$signals,
                          dim = c(1, dim(data$signals)))
    dimnames(data$signals) <- list(
      epoch = "1",
      time = sigtime,
      electrode = elecs,
      frequency = data$freq_info$freqs
    )
  }

  # Remove edges of the TFR'd data, where edge effects would be expected.
  edge_mat <- remove_edges(sigtime,
                           data$freq_info$morlet_resolution$sigma_t)

  if (keep_trials == TRUE) {
    dims <-
      which(names(dimnames(data$signals)) %in% c("time", "frequency"))
    data$signals <- sweep(data$signals,
                          dims,
                          edge_mat,
                          "*")
    data <- eeg_tfr(
      data$signals,
      srate = data$srate,
      events = data$events,
      chan_info = data$chan_info,
      reference = data$reference,
      timings = data$timings,
      freq_info = data$freq_info,
      dimensions = names(dimnames(data$signals)),
      epochs = data$epochs
    )
    if (verbose) {
      message("Returning single-trial data.")
    }
    return(data)
  }

  if (verbose) {
    message("Returning signal averaged over all trials.")
  }
  dims <-
    which(names(dimnames(data$signals)) %in% c("time", "frequency"))

  data$signals <- sweep(data$signals,
                        dims,
                        edge_mat,
                        "*")

  recording_id <- epochs(data)$recording[[1]]
  participant_id <- epochs(data)$participant_id[[1]]

  data <- eeg_tfr(
    data$signals,
    srate = data$srate,
    events = NULL,
    epochs = tibble::tibble(epoch = 1,
                            recording = recording_id,
                            participant_id = participant_id),
    chan_info = data$chan_info,
    reference = data$reference,
    timings = unique(data$timings),#data$timings[1:length(sigtime), c("epoch", "time")],
    freq_info = data$freq_info,
    dimensions = c("epoch",
                   "time",
                   "electrode",
                   "frequency")
  )
  data
}

#' Morlet wavelet
#'
#' Generate Morlet wavelet family
#'
#' @param frex Frequency range of interest
#' @param srate Sampling rate of signal
#' @param n_cycles Length of wavelet in cycles
#' @param n_freq number of frequencies to resolve
#' @param gauss_width Size of filter kernel in standard deviations
#' @keywords internal

morlet <- function(frex,
                   srate,
                   n_cycles,
                   n_freq,
                   gauss_width = 3) {

  # calculate frequency and temporal std devs
  sigma_t <- morlet_res(frex, n_cycles)$sigma_t
  max_sd_t <- max(sigma_t)
  tstep <- 1 / srate

  # round the max SD to the next biggest number divisible by tstep
  round_sd <- max_sd_t + (tstep - (max_sd_t %% tstep))

  # calculate length of kernel as 3 standard deviations
  # TO DO - change this to do it for every frequency separately
  wavtime <-
    seq(
      -round_sd * gauss_width,
      round_sd * gauss_width,
      by = tstep
    )

  t_by_f <-
    matrix(wavtime,
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
  data.frame(frequency = frex,
             sigma_f,
             sigma_t,
             n_cycles)
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
    stats::fft(signal_n)
  } else {
    if (nrow(signal) < n) {
      signal_n <- matrix(0,
                         nrow = n,
                         ncol = ncol(signal))
      signal_n[1:nrow(signal), ] <- signal
    } else {
      signal_n <- signal[1:n, ]
    }
    stats::mvfft(as.matrix(signal_n))
  }

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

#' Convert Fourier output to power, phase, or ITC as requested.
#'
#' @param data Fourier coefficients from `eeg_tfr`
#' @param output What output is desired - "power", "phase"
#' @keywords internal
convert_tfr <- function(data,
                        output) {

  if (!is.eeg_tfr(data)) {
    stop("This function requires an `eeg_tfr` object as input.")
  }

  if (is.complex(data$signals[[1]])) {
    data$signals <- switch(output,
                           "power" = abs(data$signals) ^ 2,
                           "phase" = atan2(Im(data$signals),
                                           Re(data$signals)),
                           #"itc" = calc_ITC(data),
                           "fourier" = data$signals)
  } else {
    warning("Data is not complex, returning original data.")
  }
  data
}

#' Calculate inter-trial coherence
#'
#' Calculates inter-trial coherence (ITC), a measure of phase consistency across
#' single trial data. Input data must be provided as complex Fourier
#' coefficients within an `eeg_tfr` object
#'
#' @param data An `eeg_tfr` object
#' @return An `eeg_tfr` object
#' @export
compute_itc <- function(data) {

  if (!is.eeg_tfr(data)) {
    stop("This function requires an `eeg_tfr` object as input.")
  }

  if (!is.complex(data$signals[[1]])) {
    stop("Data is not complex, returning original data.")
  }

  orig_dimnames <- dimnames(data$signals)
  orig_dim <- dim(data$signals)
  n_epochs <- orig_dim[1]

  if (n_epochs == 1) {
    stop("Only one epoch, so ITC cannot be calculated.")
  }

  data$signals <- data$signals / abs(data$signals)
  data$signals <- abs(apply(data$signals,
                            c(2, 3, 4),
                            sum))
  data$signals <- data$signals / n_epochs
  dim(data$signals) <- c(1,
                         orig_dim[2:4])
  dimnames(data$signals) <- c(epoch = list("1"),
                              orig_dimnames[2:4])
  data$timings <- data$timings[data$timings$epoch == 1, ]
  epochs(data) <- epochs(data)[epochs(data)$epoch == 1,
                               c("epoch", "participant_id")]
  data$freq_info$output <- "itc"
  class(data) <- c("tfr_average",
                   "eeg_tfr")
  data
}

#' Normalise Morlet wavelets
#'
#' Normalise the FFT'd Morlet wavelet family as recommended by Mike X Cohen,
#' dividing each wavelet by its absolute maximum. This should result in each
#' frequency being passed at unit amplitude, and the resulting convolution with
#' the signal should return units approximately on the original scale (i.e. uV^2
#' / Hz). an alternative would be Frobenius norm?
#'
#' @param mf_zp A zero-padded, FFT'd morlet wavelet family
#' @param n_freq Number of frequencies
#' @keywords internal
wavelet_norm <- function(mf_zp, n_freq) {
  mf_zp_maxes <- apply(abs(mf_zp),
                       2,
                       which.max)
  mf_zp_maxes <- lapply(seq_along(mf_zp_maxes),
                        function(x)
                          mf_zp[mf_zp_maxes[[x]], x])
  norm_mf <- lapply(seq_along(mf_zp_maxes),
                    function(x)
                      mf_zp[, x] / mf_zp_maxes[[x]])
  norm_mf <- matrix(unlist(norm_mf),
                    ncol = n_freq)
  norm_mf
}

#' @noRd
parse_frex <- function(foi,
                       n_freq,
                       verbose,
                       spacing = NULL) {

  if (length(foi) > 2) {
    stop("No more than two frequencies should be specified.")
  } else if (length(foi) == 2) {
    foi <- c(min(foi),
             max(foi))
  } else {
    foi <- c(foi, foi)
    n_freq <- 1
  }

  if (identical(spacing, "log")) {
    frex <- exp(
      seq(log(foi[1]),
          log(foi[2]),
          length.out = n_freq)
    )
  } else {
    frex <- seq(foi[1],
                foi[2],
                length.out = n_freq)
  }



  if (verbose) {
    message(
      paste("Output frequencies using", spacing, "spacing:",
            paste(round(frex, 2), collapse = " ")
      )
    )
  }
  frex
}

#' @noRd
run_tf <- function(tmp,
                   norm_mf,
                   n_conv,
                   n_kern,
                   keep_trials,
                   time_sel,
                   output) {
  tmp <- conv_to_mat(tmp)
  nhfkn <- floor(n_kern / 2) + 1
  n_chans <- dim(tmp)[3]
  n_epochs <- dim(tmp)[2]
  n_freq <- dim(norm_mf)[2]
  n_times <- length(time_sel)
  all_times <- nhfkn + (time_sel - 1)
  orig_n <- nrow(tmp)
  tmp_epo <- array(complex(1),
                   dim = c(n_conv,
                           n_epochs))
  if (keep_trials) {
    if (identical(output, "power")) {
      tfr_out <- array(0,
                       dim = c(n_epochs,
                               n_times,
                               n_chans,
                               n_freq))
      for (i in 1:n_chans) {
        tmp_epo <- fft_n(tmp[, , i],
                         n_conv)

        for (ik in 1:n_epochs) {
          tfr_out[ik, , i, ] <-
            abs(stats::mvfft(norm_mf * tmp_epo[, ik],
                             inverse = TRUE)[all_times, ] / n_conv) ^ 2
        }
      }
    } else {
      tfr_out <- array(complex(1),
                       dim = c(n_epochs,
                               n_times,
                               n_chans,
                               n_freq))
      for (i in 1:n_chans) {
        tmp_epo <- fft_n(tmp[, , i, drop = FALSE], n_conv)
        for (ik in 1:n_epochs) {
          tfr_out[ik, , i, ] <-
            (stats::mvfft(norm_mf * tmp_epo[, ik],
                         inverse = TRUE)[all_times,] / n_conv)
        }
      }
    }
  } else {
    tfr_out <- array(0,
                     dim = c(n_times, n_chans, n_freq))
    for (i in 1:n_chans) {
      tmp_epo <- fft_n(tmp[, , i, drop = FALSE],
                       n_conv)

      for (ik in 1:n_epochs) {
        tfr_out[, i, ] <-
          tfr_out[, i, ] +
          abs(stats::mvfft(norm_mf * tmp_epo[, ik], inverse = TRUE)[all_times,] / n_conv) ^ 2
      }
    }
    tfr_out <- tfr_out / n_epochs
  }
  tfr_out
}



#' Perform Hanning time-frequency analysis
#'
#' Internal function for performing Hanning wavelet transforms using convolution
#' in frequency domain
#'
#' @param data Data in `eeg_epochs` format.
#' @param foi Frequencies of interest. Scalar or character vector of the lowest
#'   and highest frequency to resolve.
#' @param n_freq Number of frequencies to be resolved.
#' @param spacing Use linear or log spacing for frequencies.
#' @param n_cycles Number of cycles at each frequency.
#' @param keep_trials Keep single trials or average over them before returning.
#' @param output Sets whether output is power, phase, or fourier coefficients.
#' @param downsample Downsampling factor (integer).
#' @param demean Remove mean before transforming.
#' @param verbose Print informative messages in console.
#' @importFrom abind abind
#' @importFrom stats nextn
#' @keywords internal

tf_hanning <- function(data,
                      foi,
                      n_freq,
                      spacing,
                      n_cycles,
                      keep_trials,
                      output,
                      downsample,
                      demean = TRUE,
                      verbose) {

  if (verbose) {
    message("Computing TFR using Hanning windows convolution")
  }

  frex <- parse_frex(foi,
                     n_freq,
                     verbose,
                     spacing = spacing)

  if (downsample > 1 && verbose) {
    message("Downsampling after transformation by a factor of ", downsample)
  }

  n_freq <- length(unique(frex))
  # if a min and max n_cycles is specified, expand out to cycles per n_freq
  if (length(n_cycles) == 2) {
    if (identical(spacing, "linear")) {
      n_cycles <- seq(n_cycles[1],
                      n_cycles[2],
                      length.out = n_freq)
    } else if (identical(spacing, "log")) {
      n_cycles <- exp(
        seq(log(n_cycles[1]),
            log(n_cycles[2]),
            length.out = n_freq)
      )
    }
  } else if (length(n_cycles) > 2) {
    stop("n_cycles should be a vector of length 1 or length 2.")
  }

  #de-mean each epoch
  if (demean) {
    data <- rm_baseline(data,
                        verbose = verbose)
  }

  elecs <- names(data$signals)
  #fft_points <- length(unique(data$timings$time))
  sigtime <- unique(data$timings$time)

  #Create a family of morlet wavelets (unscaled)
   hann_windows <-
     hann_family(
       frex = frex,
       srate = data$srate,
       n_cycles = n_cycles
       )

  # Create a list of metadata about the TFR
  data$freq_info <- list(
    freqs = frex,
    window_resolution = morlet_res(frex,
                                   n_cycles),
    morlet_res = "compatibility only, please check window_resolution",
    method = "hanning",
    output = output,
    baseline = "none"
  )

  #n_kern <- nrow(hann_windows)
  max_length <- length(unique(data$timings$time))
  n_kern <- length(hann_windows[[1]])
  n_conv <- max_length + n_kern - 1
  n_conv <- stats::nextn(n_conv, 2)

  # zero-pad and run FFTs on morlets
   norm_mf <- lapply(
     hann_windows,
     fft_n,
     n = n_conv
   )

   norm_mf <- matrix(unlist(norm_mf),
                     ncol = length(norm_mf))
                    #byrow = TRUE)


  # # Normalise wavelets for FFT (as suggested by Mike X. Cohen):
   norm_mf <- wavelet_norm(norm_mf,
                           n_freq)

  # Run the FFT convolutions on each individual trial

  # generate a vector for selection of specific timepoints
  time_sel <- seq(1,
                  length(unique(data$timings$time)),
                  by = downsample)

  data$signals <- run_tf(data,
                         norm_mf,
                         n_conv,
                         n_kern,
                         keep_trials,
                         time_sel,
                         output)

  sigtime <- sigtime[time_sel]
  data$timings <- data$timings[data$timings$time %in% sigtime,]

  if (keep_trials == TRUE) {
    dimnames(data$signals) <- list(
      epoch = unique(data$timings$epoch),
      time = sigtime,
      electrode = elecs,
      frequency = data$freq_info$freqs
    )
  } else {
    data$signals <- array(data$signals,
                          dim = c(1, dim(data$signals)))
    dimnames(data$signals) <- list(
      epoch = "1",
      time = sigtime,
      electrode = elecs,
      frequency = data$freq_info$freqs
    )
  }

  # Remove edges of the TFR'd data, where edge effects would be expected.
  edge_mat <- remove_edges(sigtime,
                           data$freq_info$window_resolution$sigma_t)

  if (keep_trials == TRUE) {
    dims <-
      which(names(dimnames(data$signals)) %in% c("time", "frequency"))
    data$signals <- sweep(data$signals,
                          dims,
                          edge_mat,
                          "*")
    data <- eeg_tfr(
      data$signals,
      srate = data$srate,
      events = data$events,
      chan_info = data$chan_info,
      reference = data$reference,
      timings = data$timings,
      freq_info = data$freq_info,
      dimensions = names(dimnames(data$signals)),
      epochs = data$epochs
    )
    if (verbose) {
      message("Returning single-trial data.")
    }
    return(data)
  }

  if (verbose) {
    message("Returning signal averaged over all trials.")
  }
  dims <-
    which(names(dimnames(data$signals)) %in% c("time", "frequency"))

  data$signals <- sweep(data$signals,
                        dims,
                        edge_mat,
                        "*")

  recording_id <- epochs(data)$recording[[1]]
  participant_id <- epochs(data)$participant_id[[1]]

  data <- eeg_tfr(
    data$signals,
    srate = data$srate,
    events = NULL,
    epochs = tibble::tibble(epoch = 1,
                            recording = recording_id,
                            participant_id = participant_id),
    chan_info = data$chan_info,
    reference = data$reference,
    timings = unique(data$timings),
    freq_info = data$freq_info,
    dimensions = c("epoch",
                   "time",
                   "electrode",
                   "frequency")
  )
  data
}


hann_family <- function(frex,
                        srate,
                        n_cycles) {

  # calculate frequency and temporal std devs
  win_sizes <- win_samples(frex = frex,
                           n_cycles = n_cycles,
                           srate = srate)

  max_win <- max(win_sizes) * 3

    windows <- lapply(
    seq_along(frex),
    function(x) {
      current_win <- win_sizes[[x]]
      window <- signal::hanning(win_sizes[[x]])
      pi_seq <- seq(from = -(current_win - 1) / 2,
                    to = (current_win - 1)/ 2,
                    length.out = current_win)

      pi_seq <- pi_seq * 2 * pi /srate
      prepad <- ceiling((max_win - length(window)) / 2)
      postpad <- floor((max_win - length(window)) / 2)
      win_times <- pi_seq * frex[[x]]
      final_win <- complex(
         real = c(rep(0, prepad),
                  window * cos(win_times),
                  rep(0, postpad)),
         imaginary = c(rep(0, prepad),
                       window * sin(win_times),
                       rep(0, postpad))
         )
      final_win
      })
  windows
}

win_samples <- function(frex,
                        n_cycles,
                        srate) {

  samps_sec <- (1 / srate)
  cycle_lengths <- 1 / frex
  win_lengths <- cycle_lengths * n_cycles # in seconds
  ceiling(win_lengths / samps_sec) # round samples up
}


#' Calculate cycles
#'
#' A helper function for calculating the appropriate min-max cycles for a fixed
#' time window/frequency resolution for use with `compute_tfr`. For some
#' analyses you may wish to keep a fixed frequency resolution across the range
#' being analysed, which requires using a fixed time window. `compute_tfr`
#' expects the minimum and maximum number of cycles to be supplied. Use this
#' function to calculate the equivalent number of cycles at each frequency.
#'
#' @param time_win Time window in seconds.
#' @param frex Frequencies of interest.
#' @export
#' @return the number of cycles for each frequency of interest
#' @examples
#' cycle_calc(.5, seq(3, 30, length.out = 10))
#' no_scale_tfr <- compute_tfr(demo_epochs, foi = c(3, 30),
#'  n_cycles = range(cycle_calc(0.5, seq(3, 30, length.out = 10))),
#'   n_freq = 10)
#' @seealso \code{\link{compute_tfr}}

cycle_calc <- function(time_win,
                       frex) {
  time_win * frex
}
