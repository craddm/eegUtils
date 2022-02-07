#' Simple absolute value thresholding
#'
#' Reject data based on a simple absolute amplitude threshold. This marks any
#' timepoint from any electrode.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_data` or `eeg_epochs`.
#' @param threshold In microvolts. If one value is supplied, it will be treated
#'   as a +- value.
#' @param reject If TRUE, remove marked data immediately, otherwise mark for
#'   inspection/rejection. Defaults to FALSE.
#' @examples
#' ar_thresh(demo_epochs, c(100))
#' @return An object of class `eeg_data` or `eeg_epochs`
#' @export

ar_thresh <- function(data,
                      threshold,
                      reject = FALSE) {
  UseMethod("ar_thresh", data)
}

#' @export
ar_thresh.default <- function(data,
                              threshold,
                              reject = FALSE) {
  stop("Not implemented for objects of class ",
       class(data))
}

#' @describeIn ar_thresh Reject data using a simple threshold.
#' @export
ar_thresh.eeg_data <- function(data,
                               threshold,
                               reject = FALSE) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- check_thresh(data,
                                 threshold)
  crossed_thresh <- rowSums(crossed_thresh) == 0

  if (reject) {
    message("Removing ",
            sum(!crossed_thresh), " (",
            round(sum(!crossed_thresh) / nrow(data$signals) * 100, 2),
            "%) timepoints.")
    data$reject$timings <- data$timings[crossed_thresh, ]
    data$timings <- data$timings[crossed_thresh, ]
    data$events <- data$events[data$events$event_time %in% data$timings$time, ]

    data$signals <- data$signals[crossed_thresh, ]

  } else {
    data$reject$timings <- data$timings[crossed_thresh, ]
  }
  data
}

#' @describeIn ar_thresh Reject data using a simple threshold.
#' @export
ar_thresh.eeg_epochs <- function(data,
                                 threshold,
                                 reject = FALSE) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- check_thresh(data, threshold)

  crossed_thresh <- rowSums(crossed_thresh) == 1
  if (!any(crossed_thresh)) {
    message("No epochs contain samples above threshold.")
  } else {
    rej_epochs <- data.frame(
      epoch = unique(data$timings$epoch[crossed_thresh]),
      reason = "threshold")
    message(paste(nrow(rej_epochs),
                  "epochs contain samples above threshold."))
    if (reject) {
      message("Removing ", nrow(rej_epochs), " epochs.")
      data <- select_epochs(data,
                            epoch_no = rej_epochs,
                            keep = FALSE)
      } else {
        data$reject$epochs <- rej_epochs
        data$reject$timings <- data$timings[crossed_thresh,]
      }
  }
  data
}

#'@noRd
check_thresh <- function(data, threshold) {
  crossed_thresh <- data$signals > max(threshold) |
    data$signals < min(threshold)

  upper_thresh <- data$signals > max(threshold)
  lower_thresh <- data$signals < min(threshold)
  total_data <- prod(dim(data$signals))

  message(sum(upper_thresh),
          " (", round(sum(upper_thresh)/total_data * 100, 2), "%) ",
          "samples above ", max(threshold) , " uV threshold.")
  message(sum(lower_thresh),
          " (", round(sum(lower_thresh)/total_data * 100, 2), "%) ",
          "samples below ", min(threshold) , " uV threshold.")
  crossed_thresh
}

#' Channel statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a `eeg_data` or `eeg_epochs` object.
#' @param ... Other parameters passed to the functions.
#' @examples
#' channel_stats(demo_epochs)
#' @return A data frame with statistics for each channel.
#' @export

channel_stats <- function(data, ...) {
  UseMethod("channel_stats", data)
}

#' @describeIn channel_stats Calculate channel statistics for `eeg_data`
#'   objects.
#' @export
channel_stats.eeg_data <- function(data,
                                   ...) {

  chan_means <- colMeans(data$signals)
  chan_sds <- apply(data$signals, 2, stats::sd)
  chan_var <- chan_sds^2
  chan_kurt <- apply(data$signals, 2, kurtosis)
  chan_range <- apply(data$signals, 2, function(x) diff(range(x)))

  data.frame(electrode = names(data$signals),
             means = chan_means,
             sds = chan_sds,
             variance = chan_var,
             kurtosis = chan_kurt,
             minmax = chan_range
             )
}

#' Epoch statistics
#'
#' Calculate various statistics for each epoch in the data
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a `eeg_data` or `eeg_epochs` object.
#' @param ... Other parameters passed to the functions.
#' @examples
#' epoch_stats(demo_epochs)
#' @export

epoch_stats <- function(data,
                        ...) {
  UseMethod("epoch_stats", data)
}

#' @describeIn epoch_stats Calculate statistics for each epoch.
#' @export
epoch_stats.eeg_epochs <- function(data,
                                   ...) {
  data$signals$epoch <- data$timings$epoch
  data <- data.table::data.table(as.data.frame(data$signals))
  epoch_vars <- data[, lapply(.SD, var), by = epoch]
  epoch_kur <- data[, lapply(.SD, kurtosis), by = epoch]
  epoch_max <- data[, lapply(.SD, max), by = epoch]
  epoch_min <- data[, lapply(.SD, min), by = epoch]
  min_max <- data[, lapply(.SD, function(x) max(x) - min(x)), by = epoch]
  stats_out <- data.table::rbindlist(list(max = epoch_max,
                                          min = epoch_min,
                                          variance = epoch_vars,
                                          kurtosis = epoch_kur,
                                          minmax = min_max),
                                     idcol = "measure")
  stats_out
}

#' Calculate kurtosis
#'
#' @param data Data to calculate kurtosis for
#' @keywords internal

kurtosis <- function(data) {
  m4 <- mean((data - mean(data)) ^ 4)
  kurt <- m4 / (stats::sd(data) ^ 4) - 3
  kurt
}

#' Remove EOG using regression
#'
#' Calculates and removes the contribution of eye movements to the EEG signal
#' using least-squares regression. Specifically, it generate regression weights
#' based on EOG channels that are used to estimate how much activity eye
#' movements are responsible for across all channels.
#'
#' @param data Data to regress - `eeg_data` or `eeg_epochs`
#' @param heog Horizontal EOG channel labels
#' @param veog Vertical EOG channel labels
#' @param bipolarize Bipolarize the EOG channels. Only works when four channels
#'   are supplied (2 HEOG and 2 VEOG).
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @return An `eeg_data` or `eeg_epochs` object with corrections
#'   applied.
#' @export

ar_eogreg <- function(data,
                      heog,
                      veog,
                      bipolarize = TRUE) {
  UseMethod("ar_eogreg", data)
}

#' @rdname ar_eogreg
#' @export
ar_eogreg.eeg_data <- function(data,
                               heog,
                               veog,
                               bipolarize = TRUE) {

  eogreg(data,
         heog,
         veog,
         bipolarize)
}

#' @rdname ar_eogreg
#' @export
ar_eogreg.eeg_epochs <- function(data,
                                 heog,
                                 veog,
                                 bipolarize = TRUE) {

  eogreg(data,
         heog,
         veog,
         bipolarize)
}

#' @noRd
eogreg <- function(data,
                   heog,
                   veog,
                   bipolarize) {

  if (bipolarize) {
    EOG <- bip_EOG(data$signals, heog, veog)
  } else {
    HEOG <- data$signals[, heog, drop = TRUE]
    VEOG <- data$signals[, veog, drop = TRUE]
    EOG <- data.frame(HEOG, VEOG)
  }

  data_chans <- channel_names(data)[!channel_names(data) %in% c(heog, veog)]
  hmz <- solve(crossprod(as.matrix(EOG)),
               crossprod(as.matrix(EOG),
                         as.matrix(data$signals[, data_chans])))
  data$signals[, data_chans] <- data$signals[, data_chans] - crossprod(t(as.matrix(EOG)), hmz)
  data
}

#' @noRd
bip_EOG <- function(data,
                    HEOG,
                    VEOG) {
  HEOG <- data[, HEOG[1]] - data[, HEOG[2]]
  VEOG <- data[, VEOG[1]] - data[, VEOG[2]]
  EOG <- data.frame(HEOG, VEOG)
  EOG
}




ar_check_rejections <- function(data) {
  if (is.null(data$reject)) {
    message("Nothing currently marked for rejection.")
  } else {
    data$reject
  }
}
