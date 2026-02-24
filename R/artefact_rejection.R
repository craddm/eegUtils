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

  data$reject$timings <- data$timings[crossed_thresh, ]
  if (reject) {
    message("Removing ",
            sum(!crossed_thresh), " (",
            round(sum(!crossed_thresh) / nrow(data$signals) * 100, 2),
            "%) timepoints.")
    data$timings <- data$timings[crossed_thresh, ]
    data$events <- data$events[data$events$event_time %in% data$timings$time, ]
    data$signals <- data$signals[crossed_thresh, ]
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
      data$reject$epochs <- rbind(data$reject$epochs, rej_epochs)
      data$reject$timings <- data$timings[crossed_thresh, ]
    }
  }
  data
}

#'@noRd
check_thresh <- function(data, threshold) {
  upper_thresh <- data$signals > max(threshold)
  lower_thresh <- data$signals < min(threshold)
  total_data <- prod(dim(data$signals))

  message(sum(upper_thresh),
          " (", round(sum(upper_thresh) / total_data * 100, 2), "%) ",
          "samples above ", max(threshold), " uV threshold.")
  message(sum(lower_thresh),
          " (", round(sum(lower_thresh) / total_data * 100, 2), "%) ",
          "samples below ", min(threshold), " uV threshold.")
  upper_thresh | lower_thresh
}

#' Channel statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An `eeg_data` or `eeg_epochs` object.
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
#' @param data An `eeg_epochs` object.
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
  dt <- data.table::data.table(data$signals)
  dt[, epoch := data$timings$epoch]

  data.table::rbindlist(list(
    max      = dt[, lapply(.SD, max), by = epoch],
    min      = dt[, lapply(.SD, min), by = epoch],
    variance = dt[, lapply(.SD, var), by = epoch],
    kurtosis = dt[, lapply(.SD, kurtosis), by = epoch],
    minmax   = dt[, lapply(.SD, function(x) max(x) - min(x)), by = epoch]
  ), idcol = "measure")
}

#' Calculate kurtosis
#'
#' @param data Data to calculate kurtosis for
#' @keywords internal

kurtosis <- function(data) {
  m <- mean(data)
  m4 <- mean((data - m) ^ 4)
  m4 / stats::var(data) ^ 2 - 3
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
    eog <- bip_eog(data$signals, heog, veog)
  } else {
    heog <- data$signals[, heog, drop = TRUE]
    veog <- data$signals[, veog, drop = TRUE]
    eog <- data.frame(heog, veog)
  }

  data_chans <- channel_names(data)[!channel_names(data) %in% c(heog, veog)]
  hmz <- solve(crossprod(as.matrix(eog)),
               crossprod(as.matrix(eog),
                         as.matrix(data$signals[, data_chans])))
  data$signals[, data_chans] <- data$signals[, data_chans] - crossprod(t(as.matrix(eog)), hmz)
  data
}

#' @noRd
bip_eog <- function(data,
                    heog,
                    veog) {
  heog <- data[, heog[1]] - data[, heog[2]]
  veog <- data[, veog[1]] - data[, veog[2]]
  eog <- data.frame(heog, veog)
  eog
}




ar_check_rejections <- function(data) {
  if (is.null(data$reject)) {
    message("Nothing currently marked for rejection.")
  } else {
    data$reject
  }
}

#' Peak-to-peak amplitude threshold
#'
#' Reject data based on the peak-to-peak amplitude (max - min) within each
#' epoch or time window. This is more sensitive to slow drifts and DC shifts
#' than an absolute threshold.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_data` or `eeg_epochs`.
#' @param threshold Peak-to-peak threshold in microvolts.
#' @param reject If TRUE, remove marked epochs immediately. Defaults to FALSE.
#' @examples
#' ar_peak2peak(demo_epochs, threshold = 100)
#' @return An object of the same class as `data`.
#' @export

ar_peak2peak <- function(data,
                         threshold,
                         reject = FALSE) {
  UseMethod("ar_peak2peak", data)
}

#' @export
ar_peak2peak.default <- function(data,
                                 threshold,
                                 reject = FALSE) {
  stop("Not implemented for objects of class ",
       class(data))
}

#' @describeIn ar_peak2peak Peak-to-peak threshold for epoched data.
#' @export
ar_peak2peak.eeg_epochs <- function(data,
                                    threshold,
                                    reject = FALSE) {

  dt <- data.table::data.table(data$signals)
  dt[, epoch := data$timings$epoch]
  p2p <- dt[, lapply(.SD, function(x) max(x) - min(x)), by = epoch]
  p2p_mat <- as.matrix(p2p[, -1, with = FALSE])
  exceeded <- rowSums(p2p_mat > threshold) > 0
  bad_epochs <- p2p$epoch[exceeded]

  message(length(bad_epochs),
          " epoch(s) exceed peak-to-peak threshold of ",
          threshold, " uV.")

  if (length(bad_epochs) == 0) {
    return(data)
  }

  rej_epochs <- data.frame(epoch = bad_epochs,
                           reason = "peak-to-peak")

  if (reject) {
    message("Removing ", length(bad_epochs), " epoch(s).")
    data <- select_epochs(data,
                          epoch_no = rej_epochs$epoch,
                          keep = FALSE)
  } else {
    data$reject$epochs <- rbind(data$reject$epochs, rej_epochs)
  }
  data
}

#' @describeIn ar_peak2peak Peak-to-peak threshold for continuous data.
#' @export
ar_peak2peak.eeg_data <- function(data,
                                  threshold,
                                  reject = FALSE) {

  p2p <- apply(data$signals, 2, function(x) max(x) - min(x))
  bad_chans <- names(p2p)[p2p > threshold]

  message(length(bad_chans),
          " channel(s) exceed peak-to-peak threshold of ",
          threshold, " uV across the recording.")

  if (reject) {
    message("Removing channel(s): ", paste(bad_chans, collapse = ", "))
    data$signals <- data$signals[, !names(data$signals) %in% bad_chans,
                                 drop = FALSE]
  } else {
    data$reject$channels <- union(data$reject$channels, bad_chans)
  }
  data
}

#' Maximum absolute gradient threshold
#'
#' Reject epochs based on the maximum absolute step between consecutive samples
#' (i.e. `max(abs(diff(x)))`). This catches electrode "pops" that are
#' invisible to amplitude or peak-to-peak checks.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_epochs`.
#' @param threshold Maximum allowable step between consecutive samples, in
#'   microvolts.
#' @param reject If TRUE, remove marked epochs immediately. Defaults to FALSE.
#' @examples
#' ar_gradient(demo_epochs, threshold = 50)
#' @return An `eeg_epochs` object.
#' @export

ar_gradient <- function(data,
                        threshold,
                        reject = FALSE) {
  UseMethod("ar_gradient", data)
}

#' @export
ar_gradient.default <- function(data,
                                threshold,
                                reject = FALSE) {
  stop("Not implemented for objects of class ",
       class(data))
}

#' @describeIn ar_gradient Gradient threshold for epoched data.
#' @export
ar_gradient.eeg_epochs <- function(data,
                                   threshold,
                                   reject = FALSE) {

  dt <- data.table::data.table(data$signals)
  dt[, epoch := data$timings$epoch]

  grad <- dt[, lapply(.SD, function(x) max(abs(diff(x)))), by = epoch]
  grad_mat <- as.matrix(grad[, -1, with = FALSE])
  exceeded <- rowSums(grad_mat > threshold) > 0
  bad_epochs <- grad$epoch[exceeded]

  message(length(bad_epochs),
          " epoch(s) exceed gradient threshold of ",
          threshold, " uV.")

  if (length(bad_epochs) == 0) {
    return(data)
  }

  rej_epochs <- data.frame(epoch = bad_epochs,
                           reason = "gradient")

  if (reject) {
    message("Removing ", length(bad_epochs), " epoch(s).")
    data <- select_epochs(data,
                          epoch_no = rej_epochs$epoch,
                          keep = FALSE)
  } else {
    data$reject$epochs <- rbind(data$reject$epochs, rej_epochs)
  }
  data
}

#' Flat channel / dead channel detection
#'
#' Flag channels (or channel×epoch combinations) where the signal variance
#' falls below a threshold, indicating a disconnected or flat electrode.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_data` or `eeg_epochs`.
#' @param threshold Minimum acceptable variance (in uV^2). Channels with
#'   variance below this value are flagged. Defaults to 0.01.
#' @param reject If TRUE, remove flagged channels or epochs immediately.
#'   Defaults to FALSE.
#' @examples
#' ar_flat(demo_epochs, threshold = 0.01)
#' @return An object of the same class as `data`.
#' @export

ar_flat <- function(data,
                    threshold = 0.01,
                    reject = FALSE) {
  UseMethod("ar_flat", data)
}

#' @export
ar_flat.default <- function(data,
                            threshold = 0.01,
                            reject = FALSE) {
  stop("Not implemented for objects of class ",
       class(data))
}

#' @describeIn ar_flat Flat channel detection for continuous data.
#' @export
ar_flat.eeg_data <- function(data,
                             threshold = 0.01,
                             reject = FALSE) {

  chan_var <- apply(data$signals, 2, stats::var)
  bad_chans <- names(chan_var)[chan_var < threshold]

  message(length(bad_chans),
          " channel(s) have variance below threshold of ",
          threshold, " uV^2.")

  if (reject && length(bad_chans) > 0) {
    message("Removing channel(s): ", paste(bad_chans, collapse = ", "))
    data$signals <- data$signals[, !names(data$signals) %in% bad_chans,
                                 drop = FALSE]
  } else {
    data$reject$channels <- union(data$reject$channels, bad_chans)
  }
  data
}

#' @describeIn ar_flat Flat channel detection for epoched data.
#' @export
ar_flat.eeg_epochs <- function(data,
                               threshold = 0.01,
                               reject = FALSE) {

  dt <- data.table::data.table(data$signals)
  dt[, epoch := data$timings$epoch]

  epoch_var <- dt[, lapply(.SD, stats::var), by = epoch]
  var_mat <- as.matrix(epoch_var[, -1, with = FALSE])
  exceeded <- rowSums(var_mat < threshold) > 0
  bad_epochs <- epoch_var$epoch[exceeded]

  message(length(bad_epochs),
          " epoch(s) have channels with variance below threshold of ",
          threshold, " uV^2.")

  if (length(bad_epochs) == 0) {
    return(data)
  }

  rej_epochs <- data.frame(epoch = bad_epochs,
                           reason = "flat")

  if (reject) {
    message("Removing ", length(bad_epochs), " epoch(s).")
    data <- select_epochs(data,
                          epoch_no = rej_epochs$epoch,
                          keep = FALSE)
  } else {
    data$reject$epochs <- rbind(data$reject$epochs, rej_epochs)
  }
  data
}
