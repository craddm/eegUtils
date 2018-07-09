#' Print Values
#'
#' Print a basic summary of the contents of an \code{eeg_data} object
#'
#' @param x \code{eeg_data} object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_data <- function(x, ...) {
  elec_names <- names(x$signals)
  n_chan <- length(elec_names)
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Electrode names\t:", elec_names, "\n")
  cat("Sampling rate\t:", x$srate, "Hz\n")
}

#' Print Values
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
  cat("Epoched EEG data\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Number of epochs\t:\t", n_epochs, "\n")
  cat("Epoch limits\t\t:\t", min(unique(x$timings$time)), "-", max(unique(x$timings$time)), "seconds\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Sampling rate\t\t:\t", x$srate, " Hz\n")
}
