#' FASTER EEG artefact rejection
#'
#' Whelan et al (2011)
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_epochs}
#' @param ... Parameters passed to FASTER

eeg_FASTER <- function(data, ...) {
  UseMethod("eeg_FASTER", data)
}


#' @describeIn eeg_FASTER Run FASTER on \code{eeg_epochs}
eeg_FASTER.eeg_epochs <- function(data, ...) {


}


#' Simple thresholding
#'
#' Reject data based on a simple absolute threshold. This goes marks any timepoint from any electrode
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data} or \code{eeg_epochs}.
#' @param threshold In microvolts. If one value is supplied, it will be treated as a +- value.
#' @param reject If TRUE, remove marked data immediately, otherwise mark for inspection/rejection. Defaults to FALSE.
#' @param ... Other arguments passed to eeg_ar_thresh
#'

eeg_ar_thresh <- function(data, threshold, reject = FALSE, ...) {
  UseMethod("eeg_ar_thresh", data)
}


#' @describeIn eeg_ar_thresh Reject data using a simple threshold.
eeg_ar_thresh.eeg_data <- function(data, threshold, reject = FALSE, ...) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- data$signals > max(threshold) | data$signals < min(threshold)

  if (reject) {
    crossed_thresh <- rowSums(crossed_thresh) == 0
    data$timings <- data$timings[crossed_thresh, ]
    data$signals <- data$signals[crossed_thresh, ]
    data$events <- data$events[data$events$event_time %in% data$timings$time, ]
    data$reference$ref_data <- data$reference$ref_data[crossed_thresh, ]
    #data <- select_times(data, data$timings[crossed_thresh, ]$time) # consider creating select_timerange vs select_timepoints
  } else {
    data$reject <- crossed_thresh
  }
  data
}


#' @describeIn eeg_ar_thresh Reject data using a simple threshold.
eeg_ar_thresh.eeg_epochs <- function(data, threshold, reject = FALSE, ...) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- data$signals > max(threshold) | data$signals < min(threshold)

  if (reject) {
    crossed_thresh <- rowSums(crossed_thresh) == 1
    rej_epochs <- unique(data$timings$epoch[crossed_thresh])
    data <- select_epochs(data, rej_epochs, keep = FALSE) # consider creating select_timerange vs select_timepoints
  } else {
    data$reject <- rej_epochs
  }
  data
}
