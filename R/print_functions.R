#' Print `eeg_data` summary
#'
#' Print a basic summary of the contents of an `eeg_data` object
#'
#' @param x `eeg_data` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_data <- function(x,
                           ...) {
  elec_names <- channel_names(x)
  n_chan <- length(elec_names)
  sig_times <- round(range(x$timings$time), 3)
  cat("EEG data\n\n")
  cat("Number of channels\t:", n_chan, "\n")
  cat("Electrode names\t\t:", elec_names, "\n")
  cat("Sampling rate\t\t:", x$srate, "Hz\n")
  cat("Reference\t\t:", x$reference$ref_chans, "\n")
  cat("Signal length:", sig_times, "seconds")
  invisible(x)
}

#' Print `eeg_epochs` summary
#'
#' Print a basic summary of the contents of an `eeg_epochs` object
#'
#' @param x `eeg_epochs` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_epochs <- function(x,
                             ...) {
  elec_names <- names(x$signals)
  n_chan <- length(elec_names)
  n_epochs <- length(unique(x$timings$epoch))
  cat("Epoched EEG data\n\n")
  cat("Number of channels\t:", n_chan, "\n")
  cat("Number of epochs\t:", n_epochs, "\n")
  cat("Epoch limits\t\t:", round(min(unique(x$timings$time)), 3),
      "-", round(max(unique(x$timings$time)), 3), "seconds\n")
  cat("Electrode names\t\t:", elec_names, "\n")
  cat("Sampling rate\t\t:", x$srate, " Hz\n")
  cat("Reference\t\t:", x$reference$ref_chans, "\n")
  invisible(x)
}

#' Print `eeg_epochs` summary
#'
#' Print a basic summary of the contents of an `eeg_epochs` object
#'
#' @param x `eeg_epochs` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_ICA <- function(x,
                          ...) {
  elec_names <- names(x$signals)
  n_chan <- length(elec_names)
  n_epochs <- length(unique(x$timings$epoch))
  cat("Epoched ICA decomposition\n\n")
  cat("Number of components\t:", n_chan, "\n")
  cat("Number of epochs\t:", n_epochs, "\n")
  cat("Epoch limits\t\t:", round(min(unique(x$timings$time)), 3),
      "-", round(max(unique(x$timings$time)), 3), "seconds\n")
  cat("Sampling rate\t\t:", x$srate, " Hz\n")
  invisible(x)
}

#' Print `eeg_tfr` summary
#'
#' Print a basic summary of the contents of an `eeg_tfr` object
#'
#' @param x `eeg_tfr` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_tfr <- function(x,
                          ...) {
  elec_names <- dimnames(x$signals)[["electrode"]]
  n_chan <- length(elec_names)
  if ("epoch" %in% names(x$timings)) {
    n_epochs <- length(unique(x$timings$epoch))
  } else {
    n_epochs <- "None, averaged."
  }
  cat("Epoched EEG TFR data\n\n")
  cat("Frequency range\t\t:\t", round(x$freq_info$freqs, 2), "\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Number of epochs\t:\t", n_epochs, "\n")
  cat("Epoch limits\t\t:\t",
      round(min(unique(x$timings$time)), 3),
      "-",
      round(max(unique(x$timings$time)), 3),
      "seconds\n")
  cat("Sampling rate\t\t:\t", x$srate, " Hz\n")
  invisible(x)
}

#' Print `eeg_evoked` summary
#'
#' Print a basic summary of the contents of an `eeg_epochs` object
#'
#' @param x `eeg_epochs` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_evoked <- function(x,
                             ...) {
  elec_names <- names(x$signals)
  n_chan <- length(elec_names)
  cat("Evoked EEG data\n\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Epoch limits\t\t:",
      round(min(unique(x$timings$time)), 3),
      "-",
      round(max(unique(x$timings$time)), 3),
      "seconds\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Sampling rate\t\t:\t", x$srate, " Hz\n")
  invisible(x)
}


#' Print `eeg_stats` summary
#'
#' Print a basic summary of the contents of an `eeg_stats` object
#'
#' @param x `eeg_stats` object to be printed
#' @param ... Further arguments passed
#' @export

print.eeg_stats <- function(x, ...) {
  elec_names <- names(x$statistic)
  n_chan <- length(elec_names)
  cat("EEG Stats\n\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Epoch limits\t\t:",
      round(min(unique(x$timings)), 3),
      "-",
      round(max(unique(x$timings)), 3),
      "seconds\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Type\t\t:\t", x$method, "\n")
  invisible(x)
}

#' Print `eeg_lm` summary
#'
#' Print a basic summary of the contents of an `eeg_lm` object
#'
#' @param x `eeg_lm` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_lm <- function(x, ...) {
  elec_names <- x$chan_info$electrode
  n_chan <- length(elec_names)
  cat("EEG Linear Model\n\n")
  cat("Formula:", paste(x$formula), "\n\n")
  cat("Number of fitted channels\t:\t", n_chan, "\n")
  cat("Channel names\t\t\t:", elec_names, "\n")
  cat("Epoch limits\t\t\t:",
      round(min(unique(x$timings$time)), 3),
      "-",
      round(max(unique(x$timings$time)), 3),
      "seconds\n")
  invisible(x)
}

#' Print `eeg_group` summary
#'
#' Print a basic summary of the contents of an `eeg_group` object
#'
#' @param x `eeg_group` object to be printed
#' @param ... Further arguments passed
#' @export
print.eeg_group <- function(x, ...) {
  cat("EEG Group Data\n\n")
  n_participants <- length(unique(epochs(x)$participant_id))
  elec_names <- x$chan_info$electrode
  n_chan <- length(elec_names)
  if (inherits(x, "eeg_evoked")) {
    cat("EEG evoked data (ERPs)\n\n")
  } else if (inherits(x, "eeg_ICA")) {
    cat("EEG ICA/SSD decompositions")
  }
  cat("Number of participants\t:\t", n_participants, "\n")
  cat("Number of channels\t:\t", n_chan, "\n")
  cat("Epoch limits\t\t:",
      round(min(unique(x$timings$time)), 3),
      "-",
      round(max(unique(x$timings$time)), 3),
      "seconds\n")
  cat("Electrode names\t\t:\t", elec_names, "\n")
  cat("Sampling rate\t\t:\t", x$srate, " Hz\n")
  invisible(x)
}
