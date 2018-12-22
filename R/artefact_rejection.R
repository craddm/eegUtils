#' FASTER EEG artefact rejection
#'
#' An implementation of the FASTER artefact rejection method for EEG by Nolan,
#' Whelan & Reilly (2010) FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods. Not yet fully implemented.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param .data An object of class \code{eeg_epochs}
#' @param ... Parameters passed to FASTER
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

eeg_FASTER <- function(.data) {
  UseMethod("eeg_FASTER", .data)
}

#' @describeIn eeg_FASTER Run FASTER on \code{eeg_epochs}
#' @export
eeg_FASTER.eeg_epochs <- function(.data, ...) {

  check_ci_str(.data$chan_info)

  channels(.data) <- validate_channels(channels(.data),
                                       channel_names(.data))

  # TODO - keep a record of which trials/channels etc are removed/interpolated
  # and allow marking for inspection rather than outright rejection.

  if (is.null(.data$reference)) {
    orig_ref <- NULL
    excluded <- NULL
  } else {
    orig_ref <- .data$reference$ref_chans
    excluded <- .data$reference$excluded
  }

  orig_chan_info <- .data$chan_info
  # Re-reference to single electrode, any should be fine, Fz is arbitrary default.
  # Note - should allow user to specify in case Fz is a known bad electrode.

  # if ("Fz" %in% names(data$signals)) {
  #   data <- reref_eeg(data, ref_chans = "Fz", exclude = excluded)
  # } else {
  #   data <- reref_eeg(data, ref_chans = names(data$signals)[14], exclude = excluded)
  # }

  orig_names <- names(.data$signals)
  # Exclude ref chan from subsequent computations (may be better to alter reref_eeg...)
  data_chans <- !(orig_names %in% .data$reference$ref_chans)

  # Step 1: channel statistics
  bad_chans <- faster_chans(.data$signals[, data_chans])
  bad_chan_n <- names(.data$signals)[bad_chans]
  message(paste("Globally bad channels:",
                paste(bad_chan_n,
                      collapse = " ")))

  if (length(bad_chan_n) > 0) {

    if (is.null(.data$chan_info)) {
      warning("no chan_info, removing chans.")
      .data <- select_elecs(.data,
                           electrode = bad_chan_n,
                           keep = FALSE)
    } else {

      #Check for any bad channels that are not in the chan_info
      check_bads <- bad_chan_n %in% .data$chan_info$electrode

      # check for any bad channels that have missing coordinates
      which_bad <- .data$chan_info$electrode %in% bad_chan_n
      missing_coords <- FALSE

      if (any(which_bad)){
        missing_coords <- apply(is.na(.data$chan_info[which_bad, ]), 1, any)
      }

      missing_bads <- bad_chan_n[!check_bads | missing_coords]

      if (length(missing_bads) > 0 ) {
        bad_chan_n <- bad_chan_n[!bad_chan_n %in% missing_bads]
        warning("Missing chan_info for bad channel(s): ",
                paste0(missing_bads,
                       collapse = " "), ". Removing channels.")
        .data <- select_elecs(.data,
                              electrode = missing_bads,
                              keep = FALSE)
      }

      if (length(bad_chan_n) > 0) {
        .data <- interp_elecs(.data,
                              bad_chan_n)
      }
    }
  }

  # Step 2: epoch statistics
  bad_epochs <- faster_epochs(.data)
  bad_epochs <- unique(.data$timings$epoch)[bad_epochs]
  message(paste("Globally bad epochs:",
                paste(bad_epochs,
                      collapse = " ")))
  .data <- select_epochs(.data,
                        epoch_no = bad_epochs,
                        keep = FALSE)

  # Step 3: ICA stats (not currently implemented)

  # Step 4: Channels in Epochs
  .data <- faster_cine(.data)

  # Step 5: Grand average step (not currently implemented, probably never will be!)

  # Return to original reference, if one existed.
  if (!is.null(orig_ref)) {
    .data <- reref_eeg(.data,
                      ref_chans = orig_ref,
                      exclude = excluded)
  }

  .data$chan_info <- orig_chan_info
  .data
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

faster_epochs <- function(.data, ...) {
  chans <- channel_names(.data)
  .data <- data.table::as.data.table(.data)
  chan_means <- .data[, lapply(.SD, mean), .SDcols = chans]
  #colMeans(data$signals)
  #data$signals <- split(data$signals,
  #                      data$timings$epoch)
  epoch_range <- .data[, lapply(.SD, function(x) max(x) - min(x)),
                          .SDcols = chans,
                          by = epoch]
  epoch_range <- epoch_range[, .(Mean = rowMeans(.SD)), by = epoch]
  epoch_range <- abs(scale(epoch_range$Mean)) > 3

  epoch_diffs <- .data[, lapply(.SD, mean),
                       .SDcols = chans,
                       by = epoch][, lapply(.SD, function(x) x - mean(x)),
                                   .SDcols = chans][ ,
                                                     .(Mean = rowMeans(.SD))]
  epoch_diffs <- abs(scale(epoch_diffs$Mean)) > 3

  epoch_vars <- .data[, lapply(.SD, var), .SDcols = chans,
                      by = epoch][, apply(.SD, 1, mean),
                                  .SDcols = chans]
  epoch_vars <- abs(scale(epoch_vars)) > 3

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

faster_cine <- function(.data, ...) {

  # get xyz coords only
  xyz_coords <- .data$chan_info[, c("electrode",
                                    "cart_x",
                                    "cart_y",
                                    "cart_z")]
  # check for rows with missing values
  missing_values <- apply(xyz_coords, 1, function(x) any(is.na(x)))
  #remove any rows with missing values
  xyz_coords <- xyz_coords[!missing_values, ]

  keep_chans <- names(.data$signals) %in% xyz_coords$electrode

  epochs <- split(.data$signals,
                  .data$timings$epoch)

  # Work out which chans are bad according to FASTER in each epoch
  bad_chans <- lapply(epochs,
                      faster_epo_stat)

  # remove channel names that are for channels we have no locations for
  bad_chans <- lapply(bad_chans,
                      function(x) x[x %in% xyz_coords$electrode])

  # remove any epochs where there were no bad channels
  bad_chans <- bad_chans[lapply(bad_chans, length) > 0]

  # Get a transfer matrix for each epoch
  bad_coords <- lapply(bad_chans,
                       function(x) interp_weights(xyz_coords,
                                                  x))

  # Work out which epochs have no transfer matrix (either all good or no bad
  # chans with locations)
  # good_epochs <- vapply(bad_coords,
  #                       is.null,
  #                       FUN.VALUE = logical(1))

  bad_coords <- bad_coords[lapply(bad_coords, length) > 0]

  # If there's nothing bad in any epoch, return the data
  if (length(bad_coords) == 0) {
    return(.data)
  }

  bad_epochs <- names(bad_coords)

  new_epochs <- lapply(bad_epochs,
                       function(x) interp_chans(epochs[[x]],
                                                bad_chans[[x]],
                                                !keep_chans,
                                                bad_coords[[x]]))

  epochs <- replace(epochs, bad_epochs, new_epochs)
  epochs <- data.table::rbindlist(epochs)
  .data$signals <- as.data.frame(epochs)
  .data
}

#' @noRd
interp_weights <- function(xyz_coords, x) {

  xyz_coords <- norm_sphere(xyz_coords)
  # rads <- sqrt(rowSums(xyz_coords[, c("cart_x", "cart_y", "cart_z")] ^ 2))
  # xyz_coords[, c("cart_x", "cart_y", "cart_z")] <-
  #   xyz_coords[, c("cart_x", "cart_y", "cart_z")] / rads

  bad_coords <- xyz_coords[xyz_coords$electrode %in% x, ]

  if (nrow(bad_coords) == 0) {
    return(NULL)
  }

  good_coords <- xyz_coords[!xyz_coords$electrode %in% x, ]

  transfer_mat <- spheric_spline(good_coords[, c("cart_x", "cart_y", "cart_z")],
                                 xyz_coords[, c("cart_x", "cart_y", "cart_z")])

  transfer_mat
}

#' Quickly calculate simple Hurst exponent for a matrix
#'
#' @param data matrix of EEG signals
#' @importFrom data.table data.table
#' @keywords internal

quick_hurst <- function(.data) {
  n <- nrow(.data)
  .data <- data.table::data.table(.data)
  dat_cumsum <- .data[, lapply(.SD, cumsum)]
  rs <- dat_cumsum[, lapply(.SD, max)] - dat_cumsum[, lapply(.SD, min)]
  rs <- rs / .data[, lapply(.SD, stats::sd)]#column_sd
  as.numeric(log(rs) / log(n))
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

  crossed_thresh <- rowSums(crossed_thresh) == 1
  rej_epochs <- unique(data$timings$epoch[crossed_thresh])
  if (reject) {
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
#' @keywords internal

epoch_stats <- function(data, ...) {
  UseMethod("epoch_stats", data)
}

#' @describeIn epoch_stats Calculate statistics for each epoch.
#' @keywords internal
epoch_stats.eeg_epochs <- function(data, ...) {
  data$signals$epoch <- data$timings$epoch
  data <- data.table::data.table(as.data.frame(data$signals))
  epoch_vars <- data[, lapply(.SD, var), by = epoch]
  epoch_kur <- data[, lapply(.SD, kurtosis), by = epoch]
  epoch_max <- data[, lapply(.SD, max), by = epoch]
  epoch_min <- data[, lapply(.SD, min), by = epoch]
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
#' @keywords internal

kurtosis <- function(data) {
  m4 <- mean((data - mean(data)) ^ 4)
  kurt <- m4 / (stats::sd(data) ^ 4) - 3
  kurt
}

#' Remove EOG using regression
#'
#' Calculates and removes the contribution of eye movements to the EEG signal using
#' least-squares regression.
#'
#' @param .data data to regress - \code{eeg_data} or \code{eeg_epochs}
#' @param heog Horizontal EOG channel labels
#' @param veog Vertical EOG channel labels
#' @param bipolarize Bipolarize the EOG channels. Only works when four channels
#'   are supplied (2 HEOG and 2 VEOG).
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export

ar_eogreg <- function(.data,
                          heog,
                          veog,
                          bipolarize = TRUE) {
  UseMethod("eeg_ar_eogreg", .data)

}

#' @rdname ar_eogreg
#' @export
ar_eogreg.eeg_data <- function(.data,
                               heog,
                               veog,
                               bipolarize = TRUE) {

  eogreg(.data,
         heog,
         veog,
         bipolarize)
}

#' @rdname ar_eogreg
#' @export
ar_eogreg.eeg_epochs <- function(.data,
                                 heog,
                                 veog,
                                 bipolarize = TRUE) {

  eogreg(.data,
         heog,
         veog,
         bipolarize)

}

#' @noRd
eogreg <- function(.data,
                   heog,
                   veog,
                   bipolarize) {

  if (bipolarize) {
    # HEOG <- .data$signals[, heog[1]] - .data$signals[, heog[2]]
    # VEOG <- .data$signals[, veog[1]] - .data$signals[, veog[2]]
    EOG <- bip_EOG(.data$signals, heog, veog)
  } else {
    HEOG <- .data$signals[, heog, drop = TRUE]
    VEOG <- .data$signals[, veog, drop = TRUE]
    EOG <- data.frame(HEOG, VEOG)
  }

  data_chans <- channel_names(.data)[!channel_names(.data) %in% c(heog, veog)]
  hmz <- solve(crossprod(as.matrix(EOG)),
               crossprod(as.matrix(EOG),
                         as.matrix(.data$signals[, data_chans])))
  .data$signals[, data_chans] <- .data$signals[, data_chans] - crossprod(t(as.matrix(EOG)), hmz)
  .data
}

#' @noRd
bip_EOG <- function(.data,
                    HEOG,
                    VEOG) {
  HEOG <- .data[, HEOG[1]] - .data[, HEOG[2]]
  VEOG <- .data[, VEOG[1]] - .data[, VEOG[2]]
  EOG <- data.frame(HEOG, VEOG)
  EOG
}
