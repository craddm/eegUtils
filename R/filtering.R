#' Filter EEG data
#'
#' Perform IIR or FIR filtering on input EEG data of class \code{eeg_data} or
#' \code{eeg_epochs}. WARNING: with epoched data, epoch boundaries are currently
#' ignored, which can result in minor edge artifacts.
#'
#' low_freq and high_freq are the low and high cutoff frequencies. Pass low freq
#' or high freq alone to perform high-pass or low-pass filtering respectively.
#' For band-pass or band-stop filters, pass both low_freq and high_freq.
#'
#' If low_freq < high_freq, bandpass filtering is performed.
#'
#' If low_freq > high_freq, bandstop filtering is performed.
#'
#' Note that it is recommended to first zero-mean the signal using either
#' channel means or by-channel epoch means.
#'
#' The function allows parallelization using the `future` package, e.g. using
#' `plan(multiprocess)`
#'
#' @section FIR versus IIR filtering: Finite Impulse Response (FIR) filtering is
#'   performed using an overlap-add FFT method. Note that this only performs a
#'   single-pass; the data is shifted back in time by the group delay of the
#'   filter to compensate for the phase delay imposed by the linear filtering
#'   process. Infinite Impulse Response (IIR) filtering is performed using a
#'   two-pass (once forwards, once reversed) method to correct for phase
#'   alignment.
#'
#' @examples
#' plot_psd(eeg_filter(demo_epochs, low_freq = 1, high_freq = 30))
#' plot_psd(eeg_filter(demo_epochs, low_freq = 12, high_freq = 8))
#' plot_psd(eeg_filter(demo_epochs, low_freq = 12, high_freq = 8, method = "iir"))
#'
#' @param .data An \code{eeg_data} or \code{eeg_epochs} object to be filtered
#' @param ... Additional parameters.
#' @export

eeg_filter <- function(.data, ...) {
  UseMethod("eeg_filter", .data)
}

#' @param low_freq Low cutoff frequency.
#' @param high_freq High cutoff frequency.
#' @param filter_order Defaults to "auto", which automatically estimates filter order for the specified filter characteristics (defaults to 4 if method = "iir").
#' @param trans_bw Transition bandwidth of the filter. "auto" or an integer.
#'   "auto" attempts to determine a suitable transition bandwidth using the
#'   heuristic given below. Ignored if method = "iir".
#' @param method "fir" (Finite Impulse Response) or "iir" (Infinite Impulse
#'   Response). Defaults to "fir".
#' @param window Windowing function to use (FIR filtering only). Defaults to
#'   "hamming"; currently only "hamming" available.
#' @param demean Remove DC component (i.e. channel/epoch mean) before filtering. Defaults to TRUE.
#' @rdname eeg_filter
#' @export
eeg_filter.eeg_data <- function(.data,
                                low_freq = NULL,
                                high_freq = NULL,
                                filter_order = "auto",
                                trans_bw = "auto",
                                method = "fir",
                                window = "hamming",
                                demean = TRUE,
                                ...) {

  filt_pars <- parse_filt_freqs(low_freq,
                                high_freq,
                                .data$srate,
                                method)

  if (identical(method, "iir")) {

    if (identical(filter_order, "auto")) {
      filter_order <- 2
    }

    message(paste("Effective filter order:",
                  filter_order * 2,
                  "(two-pass)"))

  } else if (identical(method, "fir")) {

    if (identical(trans_bw, "auto")) {
      trans_bw <- est_tbw(filt_pars,
                          .data$srate)
    }

    message(paste("Transition bandwidth:", min(trans_bw), "Hz"))

    if (identical(filter_order, "auto")) {
      filter_order <- est_filt_order(window,
                                     trans_bw,
                                     srate = .data$srate)
    }

    message(paste("Filter order:", filter_order))

  } else {
    stop("Unknown method; valid methods are 'fir' and 'iir'.")
  }

  filt_coef <- filter_coefs(method,
                            filt_pars,
                            filter_order,
                            window)

  if (demean) {
    .data <- rm_baseline(.data) # remove DC component
  }

  if (identical(method, "iir")) {
    .data <- run_iir_n(.data,
                       filt_coef)
  } else {
    .data <- run_fir(.data,
                     filt_coef,
                     filter_order)
  }
  .data$signals <- tibble::as_tibble(.data$signals)
  .data
}


#' @rdname eeg_filter
#' @export
eeg_filter.eeg_epochs <- function(.data,
                                  low_freq = NULL,
                                  high_freq = NULL,
                                  filter_order = "auto",
                                  trans_bw = "auto",
                                  method = "fir",
                                  window = "hamming",
                                  demean = TRUE,
                                  ...) {

  filt_pars <- parse_filt_freqs(low_freq,
                                high_freq,
                                .data$srate,
                                method)

  if (identical(method, "iir")) {

    if (identical(filter_order, "auto")) {
      filter_order <- 2
    }

    message(paste("Effective filter order:",
                  filter_order * 2,
                  "(two-pass)"))

  } else if (identical(method, "fir")) {

    if (identical(trans_bw, "auto")) {
      trans_bw <- est_tbw(filt_pars,
                          .data$srate)
    }

    message(paste("Transition bandwidth:", min(trans_bw), "Hz"))

    if (identical(filter_order, "auto")) {
      filter_order <- est_filt_order(window,
                                     trans_bw,
                                     srate = .data$srate)
    }

    message(paste("Filter order:", filter_order))

  } else {
    stop("Unknown method; valid methods are 'fir' and 'iir'.")
  }

  filt_coef <- filter_coefs(method,
                            filt_pars,
                            filter_order,
                            window)

  if (demean) {
    .data <- rm_baseline(.data) # remove DC component
  }
  if (identical(method, "iir")) {
    .data <- run_iir_n(.data,
                       filt_coef)
  } else {
    .data <- run_fir(.data,
                     filt_coef,
                     filter_order)
  }
  .data$signals <- tibble::as_tibble(.data$signals)
  .data
}


#' Run FIR filter using overlap-add FFT
#'
#' @param .data Data to be filtered.
#' @param filt_coef Filter coefficients
#' @param filter_order Order of filter
#' @importFrom future.apply future_lapply
#' @importFrom purrr map_df
#' @importFrom tibble as_tibble
#' @keywords internal
run_fir <- function(.data,
                    filt_coef,
                    filter_order) {

   fft_length <- length(filt_coef) * 2 - 1
   # pad the signals with zeros to help with edge effects
   .data$signals <- purrr::map_df(.data$signals,
                                  ~pad(.,
                                       fft_length))
   .data$signals <- future.apply::future_lapply(.data$signals,
                                                signal::fftfilt,
                                                b = filt_coef,
                                                n = fft_length)
   # fftfilt filters once and thus shifts everything in time by the group delay
   # of the filter (half the filter order). Here we correct for both the
   # padding and the group delay
   .data$signals <- purrr::map_df(.data$signals,
                                  ~fix_grpdelay(.,
                                                fft_length,
                                                filter_order / 2))
   .data$signals <- tibble::as_tibble(.data$signals)
   .data
}

run_iir_n <- function(.data,
                      filt_coef) {

  .data$signals <- future.apply::future_lapply(.data$signals,
                                               signal::filtfilt,
                                               filt = filt_coef,
                                               a = 1)
  .data$signals <- tibble::as_tibble(.data$signals)
  .data
}


#' Parse filter frequency input
#'
#' Parses the frequencies input by the user, converting them to a fraction of
#' the sampling rate and setting the filter type (low-pass, high-pass,
#' band-pass, band-stop) appropriately.
#'
#' @param low_freq low frequency cutoff (Hz)
#' @param high_freq High frequency cutoff (Hz)
#' @param srate Sampling rate (Hz)
#' @param method "iir" or "fir" method.
#' @keywords internal
parse_filt_freqs <- function(low_freq,
                             high_freq,
                             srate,
                             method) {

  if (length(low_freq) > 1 | length(high_freq) > 1) {
    stop("Only one number should be passed to low_freq or high_freq.")
  }

  if (is.null(low_freq) & is.null(high_freq)) {
    stop("At least one frequency must be specified.")
  }

  if (is.null(low_freq)) {
    filt_type <- "low"
    message(paste("Low-pass",
                  toupper(method),
                  "filter at",
                  high_freq, "Hz"))
    W <- high_freq / srate
  } else if (is.null(high_freq)) {
    filt_type <- "high"
    message(paste("High-pass",
                  toupper(method),
                  "filter at",
                  low_freq,
                  "Hz"))
    W <- low_freq / srate
  } else if (low_freq > high_freq) {
    filt_type <- "stop"
    message(paste("Band-stop",
                  toupper(method),
                  "filter from",
                  high_freq,
                  "-",
                  low_freq, "Hz."))
    W <- c(high_freq / srate,
           low_freq / srate)
    orig_lf <- low_freq
    low_freq <- high_freq
    high_freq <- orig_lf
  } else if (low_freq < high_freq) {
    filt_type <- "pass"
    message(paste("Band-pass",
                  toupper(method),
                  "filter from",
                  low_freq, "-",
                  high_freq, "Hz"))
    W <- c(low_freq / (srate),
           high_freq / (srate))
  }
  list("low_freq" = low_freq,
       "high_freq" = high_freq,
       "W" = W,
       "filt_type" = filt_type)
}

#' Estimate transition bandwidth
#'
#' @param filt_pars Parsed filter information from parse_filt_freqs
#' @param srate Sampling rate (Hz)
#' @keywords internal
est_tbw <- function(filt_pars,
                    srate) {

  lbw <- NA
  hbw <- NA
  max_tbw <- NA

  # stop the max transition bandwidth from being too large for notch filters
  if (identical(filt_pars$filt_type, "stop")) {
    max_tbw <- abs(filt_pars$low_freq - filt_pars$high_freq) / 2
  }

  #' Note these are the same heuristics used by MNE-python.
  if (!is.null(filt_pars$low_freq)) {
    lbw <- min(max(filt_pars$low_freq * 0.25, 2),
               filt_pars$low_freq)
  }

  if (!is.null(filt_pars$high_freq)) {
    hbw <- min(max(filt_pars$high_freq * 0.25, 2),
               srate / 2 - filt_pars$high_freq)
  }

  tbw <- min(c(lbw, hbw, max_tbw), na.rm = TRUE)
  tbw
}

#' Estimate filter order
#' @param window Window to use (string)
#' @param tbw Transition bandwidth.
#' @param srate Sampling rate (Hz)
#' @keywords internal
est_filt_order <- function(window,
                           tbw,
                           srate) {

  if (identical(window, "hamming")) {
    filt_ord <- 3.3 * srate
  } else{
    stop("Only hamming windows currently implemented.")
  }

  filt_ord <- filt_ord / min(tbw)
  filt_ord <- max(round(filt_ord, 1))
  filt_ord <- ceiling(filt_ord / 2) * 2
  filt_ord
}


#' Generate filter coefficients
#'
#' Generate filter coefficients for an IIR or FIR filter.
#'
#' @param method IIR or FIR
#' @param filt_pars output of parse_filt_freqs
#' @param filter_order order of the filter in samples
#' @param window Ignored for IIR
#' @keywords internal
filter_coefs <- function(method,
                         filt_pars,
                         filter_order,
                         window) {

  if (identical(method, "iir")) {
    coefs <- signal::butter(filter_order,
                            filt_pars$W * 2,
                            filt_pars$filt_type)
    return(coefs)
  } else {
    win <- select_window(window,
                         filter_order)

    coefs <- list()
    for (i in 1:length(filt_pars$W)) {
      coefs[[i]] <- filt_kernel(filter_order,
                                filt_pars$W[[i]],
                                win)
      if (identical(filt_pars$filt_type,
                    "high")) {
        coefs[[i]] <- spec_inv(coefs[[i]])
      }
      if (i == 2) {
        coefs[[i]] <- spec_inv(coefs[[i]])
      }
    }

    if (length(coefs) == 2) {
      coefs <- Reduce("+", coefs)
      if (identical(filt_pars$filt_type, "pass")) {
        coefs <- spec_inv(coefs)
      }
    }
  }
  unlist(coefs)
}

#' Create windowing function
#'
#' Create a windowing function for use in creating a windowed-sinc kernel
#'
#' @param type Window function to apply
#' @param m Filter order
#' @param a alpha/beta to be used for some window functions
#' @keywords internal

select_window <- function(type,
                          m,
                          a = NULL) {

  m <- m + 1
  w <- switch(type,
              "bartlett" = signal::bartlett(m),
              "hann" = signal::hanning(m),
              "hamming" = signal::hamming(m),
              "blackman" = signal::blackman(m),
              "kaiser" = signal::kaiser(m, 5.653))
  w
}


#' Create windowed-sinc filter kernel
#'
#' Supplied with the length of the filter in samples, the cutoff frequency, and
#' the windowing function, returns a windowed-sinc filter kernel for use in FIR
#' filtering.
#'
#' @param filt_order Length of the filter in samples
#' @param cut_off Cutoff frequency (normalized as a fraction of sampling rate)
#' @param w Window
#' @keywords internal

filt_kernel <- function(filt_order,
                        cut_off,
                        w) {

  m <- zero_vec(filt_order + 1)

  sinc_fun <- function(m, fc) {
    ifelse(m == 0,
           2 * pi * fc,
           sin(2 * pi * fc * m) / m)
  }

  # multiply the sinc function by the window to create a windowed sinc kernel
  filt_kern <- w * sinc_fun(m, cut_off)
  # Normalize to unity gain
  filt_kern <- filt_kern / sum(filt_kern)
  filt_kern
}

#' Spectral inversion
#'
#' Convert high-pass to low-pass, band-pass to band-stop, and vice versa.
#'
#' @param filt_kern Filter kernel to be inverted
#' @keywords internal
spec_inv <- function(filt_kern) {
  filt_kern <- -filt_kern
  midpoint <- (length(filt_kern) + 1) / 2
  filt_kern[midpoint] <- filt_kern[midpoint] + 1
  filt_kern
}

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
      data <- sweep(data,
                    2,
                    colMeans(data))
    } else {
      data <- data - mean(data)
    }

  } else if (low_freq > high_freq) {
    filt_type <- "stop"
    message(sprintf("Band-stop IIR filter from %.4g-%.4g Hz",
                    high_freq, low_freq))
    W <- c(high_freq / (srate / 2),
           low_freq / (srate / 2))

    if (length(dim(data)) > 1) {
      data <- sweep(data,
                    2,
                    colMeans(data))
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
    data <- data[,
                 lapply(.SD,
                        function(x) signal::filtfilt(filt_coef, x)),
                 by = epoch]
  } else {
    data <- data[,
                 lapply(.SD,
                        function(x) signal::filtfilt(filt_coef, x))]
  }
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
#' @keywords internal

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


