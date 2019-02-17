#' Print eeg_data summary
#'
#' Print a basic summary of the contents of an \code{eeg_data} object
#'
#' @param x \code{eeg_data} object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_data <- function(x,
                           ...) {
  elec_names <- channel_names(x)
  n_chan <- length(elec_names)
  cat("EEG data\n\n")
  cat("Number of channels\t:", n_chan, "\n")
  cat("Electrode names\t\t:", elec_names, "\n")
  cat("Sampling rate\t\t:", x$srate, "Hz\n")
  cat("Reference\t\t:", x$reference$ref_chans, "\n")
}

#' Print eeg_epochs summary
#'
#' Print a basic summary of the contents of an \code{eeg_epochs} object
#'
#' @param x \code{eeg_epochs} object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_epochs <- function(x, ...) {
  elec_names <- names(x$signals)
  n_chan <- length(elec_names)
  n_epochs <- length(unique(x$timings$epoch))
  cat("Epoched EEG data\n\n")
  cat("Number of channels\t:", n_chan, "\n")
  cat("Number of epochs\t\t:", n_epochs, "\n")
  cat("Epoch limits\t\t:", min(unique(x$timings$time)), "-", max(unique(x$timings$time)), "seconds\n")
  cat("Electrode names\t\t:", elec_names, "\n")
  cat("Sampling rate\t\t:", x$srate, " Hz\n")
  cat("Reference\t\t:", x$reference$ref_chans, "\n")
}

#' Print Values
#'
#' Print a basic summary of the contents of an \code{eeg_tfr} object
#'
#' @param x \code{eeg_tfr} object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_tfr <- function(x, ...) {
  elec_names <- dimnames(x$signals)[[2]]
  n_chan <- length(elec_names)
  if ("epoch" %in% names(x$timings)) {
    n_epochs <- length(unique(x$timings$epoch))
  } else {
    n_epochs <- "None, averaged."
  }
  cat("Epoched EEG TFR data\n\n")
  cat("Frequency range\t\t:\t", x$freq_info$freqs, "\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Number of epochs\t:\t", n_epochs, "\n")
  cat("Epoch limits\t\t:\t",
      min(unique(x$timings$time)),
      "-",
      max(unique(x$timings$time)),
      "seconds\n")
  cat("Sampling rate\t\t:\t", x$srate, " Hz\n")
}

#' Print eeg_evoked summary
#'
#' Print a basic summary of the contents of an \code{eeg_epochs} object
#'
#' @param x \code{eeg_epochs} object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_evoked <- function(x, ...) {
  elec_names <- names(x$signals)
  n_chan <- length(elec_names)
  cat("Evoked EEG data\n\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Epoch limits\t\t:\t", min(unique(x$timings$time)), "-", max(unique(x$timings$time)), "seconds\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Sampling rate\t\t:\t", x$srate, " Hz\n")
}

#' Print eeg_data summary
#'
#' Print a basic summary of the contents of an \code{eeg_data} object
#'
#' @param object \code{eeg_data} object to be printed
#' @param ... Further arguments passed
#' @export

summary.eeg_data <- function(object, ...) {
  elec_names <- names(object$signals)
  n_chan <- length(elec_names)
  cat("Epoched EEG data\n\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Electrode names\t:", elec_names, "\n")
  cat("Sampling rate\t:", object$srate, "Hz\n")
}
