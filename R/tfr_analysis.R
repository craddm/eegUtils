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
#' @param verbose Print informative messages in console.
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
                                   verbose = TRUE,
                                   ...) {


  tfr_obj <- switch(method,
                    "morlet" = tf_morlet(data = data,
                                         foi = foi,
                                         n_freq = n_freq,
                                         n_cycles = n_cycles,
                                         keep_trials = keep_trials,
                                         output = output,
                                         downsample = downsample,
                                         verbose = verbose),
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
                                   verbose = TRUE,
                                   ...) {
  switch(method,
         "morlet" = tf_morlet(data,
                              foi = foi,
                              n_freq = n_freq,
                              n_cycles = n_cycles,
                              keep_trials = keep_trials,
                              output = output,
                              downsample = downsample,
                              verbose = verbose),
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
#' @param downsample Downsampling factor (integer).
#' @param demean Remove mean before transforming.
#' @param lang defaults to R
#' @param verbose Print informative messages in console.
#' @importFrom abind abind
#' @importFrom stats nextn
#' @keywords internal

tf_morlet <- function(data,
                      foi,
                      n_freq,
                      n_cycles,
                      keep_trials,
                      output,
                      downsample,
                      demean = TRUE,
                      lang = "R",
                      verbose) {

  frex <- parse_frex(foi,
                     n_freq,
                     verbose)

  n_freq <- length(unique(frex))
  # if a min and max n_cycles is specified, expand out to cycles per n_freq
  if (length(n_cycles) == 2) {
    n_cycles <- seq(n_cycles[1],
                    n_cycles[2],
                    length.out = n_freq)
  } else if (length(n_cycles) > 2) {
    stop("n_cycles should be a vector of length 1 or length 2.")
  }

  #de-mean each epoch
  if (demean) {
    data <- rm_baseline(data, verbose = verbose)
  }
  elecs <- names(data$signals)
  #fft_points <- length(unique(data$timings$time))
  sigtime <- unique(data$timings$time)

  #Create a family of morlet wavelets (unscaled)
  morlet_family <- morlet(frex = frex,
                          srate = data$srate,
                          n_freq = n_freq,
                          n_cycles = n_cycles
  )

  # Create a list of metadata about the TFR
  data$freq_info <- list(freqs = frex,
                         morlet_resolution = morlet_res(frex,
                                                        n_cycles),
                         method = "morlet",
                         output = output,
                         baseline = "none")

  # This is a total hack to make the rest of the code behave with eeg_evoked data
  if (is.eeg_evoked(data)) {
    data$timings$epoch <- 1
  }

  n_kern <- nrow(morlet_family)
  #

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
  time_sel <- seq(1, length(unique(data$timings$time)), by = downsample)

  data$signals <- run_tf(data,
                         norm_mf,
                         n_conv,
                         n_kern,
                         keep_trials,
                         time_sel,
                         output)
  #data$signals <- data$signals[, time_sel, , , drop = FALSE]

  # data$signals <- s3tomat(data, morlet_family)
  sigtime <- sigtime[time_sel]
  data$timings <- data$timings[data$timings$time %in% sigtime, ]
  if (keep_trials) {
    dimnames(data$signals) <- list(epoch = unique(data$timings$epoch),
                                   time = sigtime,
                                   electrode = elecs,
                                   frequency = data$freq_info$freqs)
  } else {
    dimnames(data$signals) <- list(time = sigtime,
                                   electrode = elecs,
                                   frequency = data$freq_info$freqs)
  }

  # Remove edges of the TFR'd data, where edge effects would be expected.
  edge_mat <- remove_edges(sigtime,
                           data$freq_info$morlet_resolution$sigma_t)

  if (keep_trials) {

    dims <- which(names(dimnames(data$signals)) %in% c("time", "frequency"))
    data$signals <- sweep(data$signals,
                          dims,
                          edge_mat,
                          "*")
    data <- eeg_tfr(data$signals,
                    srate = data$srate,
                    events = data$events,
                    chan_info = data$chan_info,
                    reference = data$reference,
                    timings = data$timings,
                    freq_info = data$freq_info,
                    dimensions = names(dimnames(data$signals)),
                    epochs = data$epochs)
    return(data)
  }

  dims <- which(names(dimnames(data$signals)) %in% c("time", "frequency"))
  data$signals <- sweep(data$signals,
                        dims,
                        edge_mat,
                        "*")
  data <- eeg_tfr(data$signals,
                  srate = data$srate,
                  events = NULL,
                  epochs = data$epochs[1, ],
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
#' @param n_kern Kernel length
#' @importFrom stats mvfft
#' @keywords internal

conv_mor <- function(morlet_fam,
                     signal,
                     n,
                     n_kern) {

  n_chans <- ncol(signal)

  signal <- fft_n(as.matrix(signal),
                  n)

  #tf_matrix <- do_allfft(signal, morlet_fam)
  tf_matrix <- array(1i, dim = c(nrow(signal),
                                 n_chans,
                                 ncol(morlet_fam)))

  for (i in 1:ncol(morlet_fam)) {
    tf_matrix[, , i] <- mvfft(signal * morlet_fam[, i], #tf_mat2[, i, ],
                              inverse = TRUE) / n
  }
  #tf_matrix <- aperm(tf_matrix, c(1, 3, 2))

  nHfkn <- floor(n_kern / 2) + 1
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
#' @param data Fourier coefficients from eeg_tfr
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
  norm_mf
}

parse_frex <- function(foi,
                       n_freq,
                       verbose) {

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

  if (verbose){
    message("Output frequencies: ", paste(round(frex, 2), collapse = " "))
  }
  frex
}


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
            abs(mvfft(norm_mf * tmp_epo[,ik],
                      inverse = TRUE)[all_times, ] / orig_n)^2
        }
      }
    } else {
      tfr_out <- array(complex(1),
                       dim = c(n_epochs,
                               n_times,
                               n_chans,
                               n_freq))
      for (i in 1:n_chans) {
        tmp_epo <- fft_n(tmp[, , i], n_conv)
        for (ik in 1:n_epochs) {
          tfr_out[ik,,i, ] <-
            mvfft(norm_mf * tmp_epo[, ik],
                  inverse = TRUE)[all_times,] / orig_n
        }
      }
    }
  } else {
    tfr_out <- array(0, dim = c(n_times, n_chans, n_freq))
    for (i in 1:n_chans) {
      tmp_epo <- fft_n(tmp[,,i], n_conv)
      for (ik in 1:n_epochs) {
        tfr_out[ ,i, ] <-
          tfr_out[, i, ] +
          abs(mvfft(norm_mf * tmp_epo[,ik], inverse = TRUE)[all_times,] / orig_n)^2
      }
    }
    tfr_out <- tfr_out / n_epochs
  }
  tfr_out
}
