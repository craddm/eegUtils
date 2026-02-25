#' Filter EEG data
#'
#' Perform IIR or FIR filtering on input EEG data of class `eeg_data` or
#' `eeg_epochs`. For epoched data, each epoch is filtered independently to
#' prevent filter state from bleeding across trial boundaries.
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
#' `plan(multisession)`
#'
#' @section FIR versus IIR filtering: Finite Impulse Response (FIR) filtering is
#'   performed using an overlap-add FFT method. Note that this only performs a
#'   single-pass; the data is shifted back in time by the group delay of the
#'   filter to compensate for the phase delay imposed by the linear filtering
#'   process. Infinite Impulse Response (IIR) filtering is performed using a
#'   two-pass (once forwards, once reversed) method to correct for phase
#'   alignment. Note that the Butterworth filter designs used here can become
#'   numerically unstable with only a small increase in filter order. For most
#'   purposes, use FIR filters.
#'
#' @examples
#' plot_psd(eeg_filter(demo_epochs, low_freq = 1, high_freq = 30))
#' plot_psd(eeg_filter(demo_epochs, low_freq = 12, high_freq = 8))
#' plot_psd(eeg_filter(demo_epochs, low_freq = 12, high_freq = 8, method = "iir"))
#'
#' @param data An `eeg_data` or `eeg_epochs` object to be filtered.
#' @param ... Additional parameters.
#' @return An object of the original class with signals filtered according to
#'   the user's specifications
#' @export

eeg_filter <- function(data, ...) {
  UseMethod("eeg_filter", data)
}

#' @param low_freq Low cutoff frequency.
#' @param high_freq High cutoff frequency.
#' @param filter_order Defaults to "auto", which automatically estimates filter
#'   order for the specified filter characteristics (defaults to 4 if method =
#'   "iir").
#' @param trans_bw Transition bandwidth of the filter. "auto" or an integer.
#'   "auto" attempts to determine a suitable transition bandwidth using the
#'   heuristic given below. Ignored if method = "iir".
#' @param method "fir" (Finite Impulse Response) or "iir" (Infinite Impulse
#'   Response). Defaults to "fir".
#' @param window Windowing function to use (FIR filtering only). Defaults to
#'   "hamming"; currently only "hamming" available.
#' @param demean Remove DC component (i.e. channel/epoch mean) before filtering.
#'   Defaults to TRUE.
#' @rdname eeg_filter
#' @export
eeg_filter.eeg_data <- function(data,
                                low_freq = NULL,
                                high_freq = NULL,
                                filter_order = "auto",
                                trans_bw = "auto",
                                method = "fir",
                                window = "hamming",
                                demean = TRUE,
                                ...) {

  prep <- prepare_filter(data$srate,
                         low_freq,
                         high_freq,
                         filter_order,
                         trans_bw,
                         method,
                         window)

  if (demean) {
    data <- rm_baseline(data)
  }

  if (identical(method, "iir")) {
    data$signals <- iir_filter_signals(data$signals,
                                       prep$filt_coef)
  } else {
    data$signals <- fir_filter_signals(data$signals,
                                       prep$filt_coef,
                                       prep$filter_order)
  }
  data
}


#' @rdname eeg_filter
#' @export
eeg_filter.eeg_epochs <- function(data,
                                  low_freq = NULL,
                                  high_freq = NULL,
                                  filter_order = "auto",
                                  trans_bw = "auto",
                                  method = "fir",
                                  window = "hamming",
                                  demean = TRUE,
                                  ...) {

  prep <- prepare_filter(data$srate,
                         low_freq,
                         high_freq,
                         filter_order,
                         trans_bw,
                         method,
                         window)

  if (demean) {
    data <- rm_baseline(data)
  }

  epoch_factor <- data$timings$epoch
  epoch_ids <- unique(epoch_factor)
  n_per_epoch <- sum(epoch_factor == epoch_ids[1])

  if (identical(method, "fir") && prep$filter_order >= n_per_epoch) {
    warning("Filter order (",
            prep$filter_order,
            ") exceeds epoch length (",
            n_per_epoch,
            " samples). Consider filtering before epoching, ",
            "or use method = \"iir\".")
  }

  split_sigs <- split(data$signals, epoch_factor)

  filtered <- lapply(split_sigs,
                     function(epoch_df) {
                       if (identical(method, "iir")) {
                         iir_filter_signals(epoch_df,
                                            prep$filt_coef)
                       } else {
                         fir_filter_signals(epoch_df,
                                            prep$filt_coef,
                                            prep$filter_order)
                       }
                     })

  data$signals <- tibble::as_tibble(do.call(rbind, filtered))
  row.names(data$signals) <- NULL
  data
}

#' @rdname eeg_filter
#' @export
eeg_filter.eeg_group <- function(data,
                                 low_freq = NULL,
                                 high_freq = NULL,
                                 filter_order = "auto",
                                 trans_bw = "auto",
                                 method = "fir",
                                 window = "hamming",
                                 demean = TRUE,
                                 ...) {

  stop("Filtering not supported on eeg_group objects.")
}


#' Prepare filter parameters and coefficients
#'
#' Shared validation and coefficient-generation logic used by both
#' `eeg_filter.eeg_data` and `eeg_filter.eeg_epochs`.
#'
#' @param srate Sampling rate (Hz)
#' @param low_freq Low cutoff frequency
#' @param high_freq High cutoff frequency
#' @param filter_order Filter order or "auto"
#' @param trans_bw Transition bandwidth or "auto"
#' @param method "fir" or "iir"
#' @param window Windowing function name
#' @return A list with `filt_pars`, `filter_order`, and `filt_coef`.
#' @keywords internal
prepare_filter <- function(srate,
                           low_freq,
                           high_freq,
                           filter_order,
                           trans_bw,
                           method,
                           window) {

  filt_pars <- parse_filt_freqs(low_freq,
                                high_freq,
                                srate,
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
                          srate)
    }

    message(paste("Transition bandwidth:", min(trans_bw), "Hz"))

    if (identical(filter_order, "auto")) {
      filter_order <- est_filt_order(window,
                                     trans_bw,
                                     srate = srate)
    }

    message(paste("Filter order:", filter_order))

  } else {
    stop("Unknown method; valid methods are 'fir' and 'iir'.")
  }

  filt_coef <- filter_coefs(method,
                            filt_pars,
                            filter_order,
                            window)

  list(filt_pars = filt_pars,
       filter_order = filter_order,
       filt_coef = filt_coef)
}


#' FIR filter a data.frame of signals
#'
#' Applies an FIR filter to each column using reflected-edge padding and
#' overlap-add FFT convolution. Corrects for the linear phase (group delay) of
#' the single-pass FIR filter.
#'
#' @param signals A data.frame where each column is a channel.
#' @param filt_coef Numeric vector of FIR filter coefficients.
#' @param filter_order Filter order in samples.
#' @return A tibble with filtered signals.
#' @importFrom future.apply future_lapply
#' @keywords internal
fir_filter_signals <- function(signals,
                               filt_coef,
                               filter_order) {

  n_sig <- nrow(signals)
  pad_n <- min(filter_order, n_sig - 1)
  grp_delay <- as.integer(filter_order / 2)
  # Ensure group delay correction fits within the padding
  grp_delay <- min(grp_delay, pad_n)
  fft_length <- stats::nextn(length(filt_coef) * 2 - 1)

  filtered <- future.apply::future_lapply(signals,
                                          function(x) {
                                            x_padded <- reflect_pad(x, pad_n)
                                            x_filt <- signal::fftfilt(filt_coef,
                                                                      x_padded,
                                                                      n = fft_length)
                                            start <- pad_n + 1 + grp_delay
                                            end <- pad_n + n_sig + grp_delay
                                            x_filt[start:end]
                                          },
                                          future.packages = "eegUtils")
  tibble::as_tibble(filtered)
}


#' IIR filter a data.frame of signals
#'
#' Applies a two-pass (forward-reverse) IIR filter to each column with
#' reflected-edge padding to reduce edge transients. The padding length follows
#' SciPy's convention of `3 * max(length(b), length(a))`.
#'
#' @param signals A data.frame where each column is a channel.
#' @param filt_coef An Arma object (from [signal::butter()]).
#' @return A tibble with filtered signals.
#' @importFrom future.apply future_lapply
#' @keywords internal
iir_filter_signals <- function(signals,
                               filt_coef) {

  pad_n <- 3 * max(length(filt_coef$b), length(filt_coef$a))

  filtered <- future.apply::future_lapply(signals,
                                          function(x) {
                                            n <- length(x)
                                            pn <- min(pad_n, n - 1)
                                            x_padded <- reflect_pad(x, pn)
                                            # Forward pass
                                            fwd <- as.numeric(signal::filter(filt_coef,
                                                                             x_padded))
                                            # Reverse pass
                                            rev_filt <- as.numeric(
                                              rev(signal::filter(filt_coef,
                                                                 rev(fwd)))
                                            )
                                            # Trim padding
                                            rev_filt[(pn + 1):(pn + n)]
                                          },
                                          future.packages = "eegUtils")
  tibble::as_tibble(filtered)
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

  if (length(low_freq) > 1 || length(high_freq) > 1) {
    stop("Only one number should be passed to low_freq or high_freq.")
  }

  if (is.null(low_freq) && is.null(high_freq)) {
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
  } else {
    stop("Only hamming windows currently implemented.")
  }

  filt_ord <- filt_ord / min(tbw)
  filt_ord <- max(round(filt_ord, 1))
  filt_ord <- ceiling(filt_ord / 2) * 2
  filt_ord
}


#' Generate filter coefficients
#'
#' Generate filter coefficients for an IIR or FIR filter. FIR coefficients are
#' generated using [signal::fir1()], which places the -6 dB cutoff precisely at
#' the specified frequency.
#'
#' @param method IIR or FIR
#' @param filt_pars output of parse_filt_freqs
#' @param filter_order order of the filter in samples
#' @param window Window function name (FIR only)
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
  }

  win <- switch(window,
                "hamming" = signal::hamming(filter_order + 1),
                "hann" = signal::hanning(filter_order + 1),
                "blackman" = signal::blackman(filter_order + 1),
                stop("Unsupported window type: ", window))

  coefs <- signal::fir1(filter_order,
                        filt_pars$W * 2,
                        type = filt_pars$filt_type,
                        window = win)
  as.numeric(coefs)
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
  fx <- exp(-.5 * (x / s)^2)
  fx <- fx / max(fx)
  filt_sig <- apply(data, 2,
                    function(x) {
                                 2 * Re(stats::fft(stats::fft(x) / srate * fx,
                                                   inverse = TRUE))})
  as.data.frame(filt_sig)
}
