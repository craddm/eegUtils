#' Object creator for eeg_tfr objects.
#'
#' @param data TFR transformed data
#' @param srate Sampling rate in Hz.
#' @param events Event tables
#' @param chan_info Standard channel information.
#' @param reference Reference information
#' @param timings Timing information.
#' @param freq_info Frequencies and other useful information
#' @param dimensions List of which dimension is which
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @keywords internal
eeg_tfr <- function(data,
                    srate,
                    events,
                    chan_info = NULL,
                    reference,
                    timings = NULL,
                    freq_info,
                    dimensions) {

  value <- list(signals = data,
                srate = srate,
                events = events,
                chan_info = chan_info,
                reference = reference,
                timings = timings,
                freq_info = freq_info,
                dimensions = dimensions)
  class(value) <- "eeg_tfr"
  value
}

#' Check if object is of class eeg_tfr
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object to check.
#'
is.eeg_tfr <- function(x) inherits(x, "eeg_tfr")


#' Function to create an object of class eeg_psd
#'
#' @param data PSD transformed data
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param freqs vector of frequencies
#' @param dimensions List of which dimension is which
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @noRd
eeg_psd <- function(data,
                    srate,
                    events,
                    chan_info = NULL,
                    timings = NULL,
                    freqs,
                    dimensions) {

  value <- list(signals = data,
                srate = srate,
                chan_info = chan_info,
                timings = timings,
                freqs = freqs
                )
  class(value) <- "eeg_tfr"
  value
}


#' Function to create an object of class \code{eeg_GA}
#'
#' @noRd

eeg_GA <- function(data,
                   srate,
                   chan_info,
                   timings,
                   indivs) {

  value <- list(signals = data,
                srate = srate,
                chan_info = chan_info,
                timings = timings,
                indivs = indivs)
  class(value) <- "eeg_GA"
  value
}



#' Function to create an S3 object of class \code{eeg_data}
#'
#' @noRd
validate_eeg_data <- function(data,
                         srate,
                         events = NULL,
                         chan_info = NULL,
                         timings = NULL,
                         continuous = NULL,
                         reference = NULL) {
  stopifnot(is.double(srate))
  stopifnot(is.data.frame(events))
  structure(signals = data,
            srate = srate,
            events = events)

}

#' Function to create an S3 object of class "eeg_data".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Raw data - signals from electrodes/channels.
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param continuous Whether the data is continuous or epoched.
#' @param reference Reference channel information, including names of reference
#'   channels, excluded channels etc.
#' @param events Event table
#' @export

eeg_data <- function(data,
                     srate,
                     events = NULL,
                     chan_info = NULL,
                     timings = NULL,
                     continuous = NULL,
                     reference = NULL) {
  if (srate < 1) {
    stop("Sampling rate must be above 0")
  }
  value <- list(signals = data,
                srate = srate,
                events = events,
                chan_info = chan_info,
                timings = timings,
                continuous = continuous,
                reference = reference
  )
  class(value) <- "eeg_data"
  value
}

#' Function to create an S3 object of class "eeg_evoked"
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data evoked data
#' @param chan_info Electrode locations etc
#' @param timings vector of timepoints
#' @param ... Other parameters
#' @keywords internal

eeg_evoked <- function(data,
                       chan_info,
                       timings, ...) {
  value <- list(signals = data,
                chan_info = chan_info,
                timings = timings)
  class(value) <- c("eeg_evoked")
}

#' Function to create an S3 object of class "eeg_stats".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param statistic Calculated statistic (e.g. t-statistic)
#' @param pvals calculated p-values for that statistic
#' @param chan_info String of character names for electrodes.
#' @param timings Unique timepoints remaining in the data.
#' @keywords internal

eeg_stats <- function(statistic,
                      chan_info,
                      pvals,
                      timings) {

  value <- list(statistic = statistic,
                pvals = pvals,
                chan_info = chan_info,
                timings = timings)
  class(value) <- "eeg_stats"
  value
}


eeg_ICA <- function(mixing_matrix,
                    unmixing_matrix,
                    signals,
                    timings,
                    events,
                    chan_info,
                    srate,
                    continuous) {

  value <- list(mixing_matrix = mixing_matrix,
                unmixing_matrix = unmixing_matrix,
                signals = signals,
                timings = timings,
                events = events,
                chan_info = chan_info,
                srate = srate,
                continuous = continuous)
  class(value) <- c("eeg_ICA", "eeg_epochs")
  value

}

#' Check if object is of class "eeg_data".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object to check.
#' @keywords internal

is.eeg_data <- function(x) inherits(x, "eeg_data")

#' Check if object is of class "eeg_epochs".
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object to check.
#' @keywords internal

is.eeg_epochs <- function(x) inherits(x, "eeg_epochs")

#' Check if object is of class \code{eeg_evoked}
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object to check.
#' @keywords internal
is.eeg_evoked <- function(x) inherits(x, "eeg_evoked")

#' Check if object is of class \code{eeg_stats}
#'
#' @param x Object to check.
#' @keywords internal
is.eeg_stats <- function(x) inherits(x, "eeg_stats")

#' Check if object is of class \code{eeg_ICA}
#'
#' @param x Object to check.
#' @keywords internal

is.eeg_ICA <- function(x) inherits(x, "eeg_ICA")
