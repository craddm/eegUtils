#' FASTER EEG artefact rejection
#'
#' An implementation of the FASTER artefact rejection method for EEG by Nolan,
#' Whelan & Reilly (2010) FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods. Not yet fully implemented.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_epochs}
#' @param ... Parameters passed to FASTER
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

eeg_FASTER <- function(data, ...) {
  UseMethod("eeg_FASTER", data)
}

#' @describeIn eeg_FASTER Run FASTER on \code{eeg_epochs}
#' @export
eeg_FASTER.eeg_epochs <- function(data, ...) {

  # TODO - keep a record of which trials/channels etc are removed/interpolated
  # and allow marking for inspection rather than outright rejection.

  if (is.null(data$reference)) {
    orig_ref <- NULL
    excluded <- NULL
  } else {
    orig_ref <- data$reference$ref_chans
    excluded <- data$reference$excluded
  }

  orig_chan_info <- data$chan_info
  # Re-reference to single electrode, any should be fine, Fz is arbitrary default.
  # Note - should allow user to specify in case Fz is a known bad electrode.

  # if ("Fz" %in% names(data$signals)) {
  #   data <- reref_eeg(data, ref_chans = "Fz", exclude = excluded)
  # } else {
  #   data <- reref_eeg(data, ref_chans = names(data$signals)[14], exclude = excluded)
  # }

  orig_names <- names(data$signals)
  # Exclude ref chan from subsequent computations (may be better to alter reref_eeg...)
  data_chans <- !(orig_names %in% data$reference$ref_chans)

  # Step 1: channel statistics
  bad_chans <- faster_chans(data$signals[, data_chans])
  bad_chan_n <- names(data$signals)[bad_chans]
  message(paste("Globally bad channels:",
                paste(bad_chan_n,
                      collapse = " ")))

  if (length(bad_chan_n) > 0) {

    #check_chans <- data$chan_info$electrode %in% bad_chan_n
    if (is.null(data$chan_info)) {
      warning("no chan_info, removing chans.")
      data <- select_elecs(data,
                           electrode = bad_chan_n,
                           keep = FALSE)
    } else {

      #Check for any bad channels that are not in the chan_info
      check_bads <- bad_chan_n %in% data$chan_info$electrode

      # check for any bad channels that have missing coordinates
      which_bad <- data$chan_info$electrode %in% bad_chan_n
      missing_coords <- FALSE

      if (any(which_bad)){
        missing_coords <- apply(is.na(data$chan_info[which_bad, ]), 1, any)
      }

      missing_bads <- bad_chan_n[!check_bads | missing_coords]

      if (length(missing_bads) > 0 ) {
        bad_chan_n <- bad_chan_n[!bad_chan_n %in% missing_bads]
        warning("Missing chan_info for bad channel(s): ",
                paste0(missing_bads,
                       collapse = " "), ". Removing channels.")
        data <- select_elecs(data,
                             electrode = missing_bads,
                             keep = FALSE)
      }
      data <- interp_elecs(data,
                           bad_chan_n)
    }
  }

  # Step 2: epoch statistics
  bad_epochs <- faster_epochs(data)
  bad_epochs <- unique(data$timings$epoch)[bad_epochs]
  message(paste("Globally bad epochs:",
                paste(bad_epochs,
                      collapse = " ")))
  data <- select_epochs(data,
                        epoch_no = bad_epochs,
                        keep = FALSE)

  # Step 3: ICA stats (not currently implemented)

  # Step 4: Channels in Epochs
  data <- faster_cine(data)

  # Step 5: Grand average step (not currently implemented, probably never will be!)

  # Return to original reference, if one existed.
  if (!is.null(orig_ref)) {
    data <- reref_eeg(data,
                      ref_chans = orig_ref,
                      exclude = excluded)
  }

  data$chan_info <- orig_chan_info
  data
}

#' Perform global bad channel detection for FASTER
#'
#' @param data A matrix of EEG data signals
#' @param sds Standard deviation thresholds
#' @param ... Further parameters (tbd)
#' @keywords internal

faster_chans <- function(data, sds = 3, ...) {
  chan_hurst <- scale(quick_hurst(data))
  chan_vars <- scale(apply(data,
                           2,
                           stats::var))
  chan_corrs <- scale(colMeans(abs(stats::cor(data))))
  bad_chans <- matrix(c(abs(chan_hurst) > sds,
                        abs(chan_vars) > sds,
                        abs(chan_corrs) > sds),
                      nrow = 3,
                      byrow = TRUE)
  bad_chans <- apply(bad_chans,
                     2,
                     any)
  bad_chans
}

#' Perform global bad epoch detection for FASTER
#'
#' @param data \code{eeg_epochs} object
#' @param ... Further parameters (tbd)
#' @keywords internal

faster_epochs <- function(data, ...) {
  chan_means <- colMeans(data$signals)
  data$signals <- split(data$signals,
                        data$timings$epoch)
  epoch_range <- lapply(data$signals,
                        function(x) diff(apply(x,
                                               2,
                                               range)))
  epoch_range <- rowMeans(do.call(rbind,
                                  epoch_range))
  epoch_range <- abs(scale(epoch_range)) > 3
  # epoch_range <- abs(scale(do.call("rbind",
  #                                  epoch_range))) > 3
  epoch_diffs <- lapply(data$signals,
                        function(x) {
                          apply(abs(sweep(x,
                                          2,
                                          chan_means)),
                                2,
                                mean)
                          })
  epoch_diffs <- rowMeans(do.call(rbind,
                                  epoch_diffs))
  epoch_diffs <- abs(scale(epoch_diffs)) > 3
  # epoch_diffs <- abs(scale(do.call("rbind",
                                   # epoch_diffs))) > 3
  epoch_vars <- lapply(data$signals,
                       function(x) apply(x,
                                         2,
                                         stats::var))
  epoch_vars <- rowMeans(do.call(rbind,
                                 epoch_vars))
  epoch_vars <- abs(scale(epoch_vars)) >3
  #epoch_vars <- abs(scale(do.call("rbind",
   #                               epoch_vars))) > 3
  bad_epochs <- matrix(c(rowSums(epoch_vars) > 0,
                         rowSums(epoch_range) > 0,
                         rowSums(epoch_diffs) > 0),
                       ncol = 3)
  bad_epochs <- apply(bad_epochs, 1, any)
  bad_epochs
}

#' FASTER detection of bad channels in single epochs
#'
#' @param data \code{eeg_epochs} object.
#' @param ... further parameters (tbd)
#' @keywords internal

faster_cine <- function(data, ...) {

  epochs <- split(data$signals,
                  data$timings$epoch)
  bad_chans <- lapply(epochs,
                      faster_epo_stat)
  epoch_nos <- unique(data$timings$epoch)

  good_epochs <- vapply(bad_chans,
                       purrr::is_empty,
                       FUN.VALUE = logical(1))

  # If there's nothing bad in any epoch, return the data
  if (all(good_epochs)) {
    return(data)
  }

  # go through each epoch and interpolate bad channels
  bad_epoch_nos <- epoch_nos[!good_epochs]
  bad_epochs <- lapply(bad_epoch_nos,
                       function(x) select_epochs(data,
                                                 epoch_no = x))
  bad_epochs <- lapply(seq_along(bad_epochs),
                       function(x) interp_elecs(bad_epochs[[x]],
                                                bad_chans[[x]]))

  if (length(bad_epochs) > 1) {
    bad_epochs <- do.call("eeg_combine", bad_epochs)
  } else {
    bad_epochs <- unlist(bad_epochs)
  }

  # If there are any epochs that are ok, select them, then combine with bad
  # epochs and return
  if (any(good_epochs)) {
    clean_epochs <- select_epochs(data,
                                  epoch_no = epoch_nos[good_epochs])
    data <- eeg_combine(clean_epochs,
                        bad_epochs)
    return(data)
  }

  # only get here if there were no clean epochs
  data <- bad_epochs
  data
}


#' Quickly calculate simple Hurst exponent for a matrix
#'
#' @param data matrix of EEG signals
#' @importFrom matrixStats colMaxs colMins colSds
#' @keywords internal

quick_hurst <- function(data) {
  n <- nrow(data)
  #y <- sweep(data, 2, Matrix::colMeans(data))
  s <- apply(data, 2, cumsum)
  rs <- (matrixStats::colMaxs(s) - matrixStats::colMins(s))
  rs <- rs / matrixStats::colSds(as.matrix(data))
  log(rs) / log(n)
}

#' Calculate statistics for each channel in an epoch and identify bad channels
#'
#' @param data a matrix of signals from a single epoch
#' @keywords internal

faster_epo_stat <- function(data, chan_means) {

  measures <- data.frame(vars = matrixStats::colVars(as.matrix(data)),
                         medgrad = matrixStats::colMedians(diff(as.matrix(data))),
                         range_diff = t(diff(t(matrixStats::colRanges(as.matrix(data)))))
                         #dev = sweep(data, 2, chan_means)
                         )
  # Check if any measure is above 3 standard deviations
  bad_chans <- rowSums(scale(measures) > 3) > 0
  bad_chans <- names(data)[bad_chans]
  bad_chans
}

#' Simple absolute value thresholding
#'
#' Reject data based on a simple absolute threshold. This marks any
#' timepoint from any electrode.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data} or \code{eeg_epochs}.
#' @param threshold In microvolts. If one value is supplied, it will be treated
#'   as a +- value.
#' @param reject If TRUE, remove marked data immediately, otherwise mark for
#'   inspection/rejection. Defaults to FALSE.
#' @param ... Other arguments passed to eeg_ar_thresh
#' @export

eeg_ar_thresh <- function(data,
                          threshold,
                          reject = FALSE,
                          ...) {
  UseMethod("eeg_ar_thresh", data)
}

#' @describeIn eeg_ar_thresh Reject data using a simple threshold.
#' @export
eeg_ar_thresh.eeg_data <- function(data,
                                   threshold,
                                   reject = FALSE,
                                   ...) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- data$signals > max(threshold) |
    data$signals < min(threshold)

  if (reject) {
    crossed_thresh <- rowSums(crossed_thresh) == 0
    data$timings <- data$timings[crossed_thresh, ]
    data$signals <- data$signals[crossed_thresh, ]
    data$events <- data$events[data$events$event_time %in% data$timings$time, ]
    data$reference$ref_data <- data$reference$ref_data[crossed_thresh, ]
  } else {
    data$reject <- crossed_thresh
  }
  data
}

#' @describeIn eeg_ar_thresh Reject data using a simple threshold.
#' @export
eeg_ar_thresh.eeg_epochs <- function(data, threshold, reject = FALSE, ...) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- data$signals > max(threshold) |
    data$signals < min(threshold)

  if (reject) {
    crossed_thresh <- rowSums(crossed_thresh) == 1
    rej_epochs <- unique(data$timings$epoch[crossed_thresh])
    data <- select_epochs(data, rej_epochs, keep = FALSE)
    # consider creating select_timerange vs select_timepoints
  } else {
    data$reject <- rej_epochs
  }
  data
}

#' Channel statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param ... Other parameters passed to the functions.
#' @keywords internal

channel_stats <- function(data, ...) {
  UseMethod("channel_stats", data)
}

#' @describeIn channel_stats Calculate channel statistics for \code{eeg_data}
#'   objects.
#' @keywords internal
channel_stats.eeg_data <- function(data, ...) {

  chan_means <- colMeans(data$signals)
  #chan_means <- as.data.table(data$signals)
  #chan_means <- chan_means[, lapply(.SD, mean)]
  chan_sds <- apply(data$signals, 2, stats::sd)
  chan_var <- apply(data$signals, 2, stats::var)
  chan_kurt <- apply(data$signals, 2, kurtosis)

  data.frame(electrode = names(data$signals),
             means = chan_means,
             variance = chan_var,
             kurtosis = chan_kurt
             )
}

#' Epoch statistics
#'
#' Calculate various statistics for each epoch in the data
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param ... Other parameters passed to the functions.
#' @noRd

epoch_stats <- function(data, ...) {
  UseMethod("epoch_stats", data)
}

#' @describeIn epoch_stats Calculate statistics for each epoch.
#' @noRd
epoch_stats.eeg_epochs <- function(data, ...) {
  data$signals$epoch <- data$timings$epoch
  data <- data.table::data.table(as.data.frame(data$signals))
  epoch_vars <- data[, lapply(.SD, var), by = epoch]
  epoch_kur <- data[, lapply(.SD, kurtosis), by = epoch]
  epoch_max <- data[, lapply(.SD, max), by = epoch]
  epoch_min <- data[, lapply(.SD, min), by = epoch]
  #epoch_vars <- data[, value := matrixStats::rowVars(as.matrix(.SD), na.rm = TRUE), by = epoch][, c("epoch", "value")]
  #epoch_kur <- data[, value := kurtosis(x), by = epoch][, c("epoch", "value")]
  stats_out <- data.table::rbindlist(list(max = epoch_max,
                                          min = epoch_min,
                                          variance = epoch_vars,
                                          kurtosis = epoch_kur),
                                     idcol = "measure")
  stats_out
}

#' Calculate kurtosis
#'
#' @param data Data to calculate kurtosis for
#' @noRd

kurtosis <- function(data) {
  m4 <- mean((data - mean(data)) ^ 4)
  kurt <- m4 / (stats::sd(data) ^ 4) - 3
  kurt
}

#' Regress EOG
#'
#' @param data data to regress
#' @param heog Horizontal EOG channels
#' @param veog Vertical EOG channels
#' @param bipolarize Bipolarize the EOG channels. Only works when four channels
#'   are supplied (2 HEOG and 2 VEOG).
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @noRd

eeg_ar_eogreg <- function(data, heog, veog, bipolarize = TRUE, ...) {
  UseMethod("eeg_ar_eogreg", data)

}

eeg_ar_eogreg.eeg_data <- function(data, heog, veog, bipolarize = TRUE, ...) {

  heog_only <- select_elecs(data, electrode = heog)
  veog_only <- select_elecs(data, electrode = veog)

  EOG <- data.frame(heog = NA, veog = NA)
  if (bipolarize) {
    EOG$heog <- heog_only$signals[, 1] - heog_only$signals[, 2]
    EOG$veog <- veog_only$signals[, 1] - veog_only$signals[, 2]
  } else {
    EOG$heog <- heog_only$signals[, 1]
    EOG$veog <- veog_only$signals[, 1]
  }

  hmz <- solve(EOG * t(EOG), EOG * t(data$signals))
  hmz
}
