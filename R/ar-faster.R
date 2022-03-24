#' FASTER EEG artefact rejection
#'
#' An implementation of the FASTER artefact rejection method for EEG by Nolan,
#' Whelan & Reilly (2010) FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_epochs`
#' @param ... Parameters passed to FASTER
#' @examples
#' ar_FASTER(demo_epochs)
#' @return An `eeg_epochs` object with artefact correction applied.
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

ar_FASTER <- function(data, ...) {
  UseMethod("ar_FASTER", data)
}

#' @param exclude Channels to be ignored by FASTER.
#' @param test_chans Logical. Run tests of global channel statistics
#' @param test_epochs Logical. Run tests of globally bad epochs.
#' @param test_cine Logical. Run tests for locally bad channels within epochs.
#' @describeIn ar_FASTER Run FASTER on `eeg_epochs`
#' @export
ar_FASTER.eeg_epochs <- function(data,
                                 exclude = NULL,
                                 test_chans = TRUE,
                                 test_epochs = TRUE,
                                 test_cine = TRUE,
                                 ...) {

  check_ci_str(data$chan_info)

  if (is.null(channels(data))) {
    stop("No channel information found, ar_FASTER() requires (at least some) channel locations.")
  } else {
    channels(data) <- validate_channels(channels(data),
                                        channel_names(data))
  }

  # TODO - keep a record of which trials/channels etc are removed/interpolated
  # and allow marking for inspection rather than outright rejection.

  if (is.null(data$reference)) {
    orig_ref <- NULL
    excluded <- NULL
  } else {
    orig_ref <- data$reference$ref_chans
    excluded <- data$reference$excluded
  }

  orig_chan_info <- channels(data)
  # Re-reference to single electrode, any should be fine, Fz is arbitrary default.
  # Note - should allow user to specify in case Fz is a known bad electrode.

  # if ("Fz" %in% names(data$signals)) {
  #   data <- eeg_reference(data, ref_chans = "Fz", exclude = excluded)
  # } else {
  #   data <- eeg_reference(data, ref_chans = names(data$signals)[14], exclude = excluded)
  # }

  orig_names <- channel_names(data)
  # Exclude ref chan from subsequent computations (may be better to alter
  # reref_eeg...)
  data_chans <- orig_names[!(orig_names %in% data$reference$ref_chans)]
  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      exclude <- orig_names[exclude]
    }
    message("Excluding channel(s):",
            paste(exclude, ""))
    data_chans <- data_chans[!(data_chans %in% exclude)]
  }

  # Step 1: channel statistics
  if (test_chans) {
    bad_chans <- faster_chans(data$signals[, data_chans])
    bad_chan_n <- data_chans[bad_chans]
    message(paste("Globally bad channels:",
                  paste(bad_chan_n,
                        collapse = " ")))

    if (length(bad_chan_n) > 0) {

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

        if (any(which_bad)) {
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

        if (length(bad_chan_n) > 0) {
          data <- interp_elecs(data,
                               bad_chan_n)
        }
      }
    }
  } else {
    message("Skipping global bad channel detection...")
  }

  # Step 2: epoch statistics
  if (test_epochs) {
    bad_epochs <- faster_epochs(data)
    bad_epochs <- unique(data$timings$epoch)[bad_epochs]
    message(paste("Globally bad epochs:",
                  paste(bad_epochs,
                        collapse = " ")))
    data$reject$bad_epochs <- bad_epochs
    data <- select_epochs(data,
                          epoch_no = bad_epochs,
                          keep = FALSE)
  } else {
    message("Skipping bad epoch detection")
  }
  # Step 3: ICA stats (not currently implemented)

  # Step 4: Channels in Epochs
  data <- faster_cine(data,
                      exclude)

  # Step 5: Grand average step (not currently implemented, probably never will be!)

  # Return to original reference, if one existed.
  if (!is.null(orig_ref)) {
    data <- eeg_reference(data,
                          ref_chans = orig_ref,
                          exclude = excluded,
                          verbose = FALSE)
  }

  #data$chan_info <- orig_chan_info
  data
}

#' Perform global bad channel detection for FASTER
#'
#' @param data A matrix of EEG data signals
#' @param sds Standard deviation thresholds
#' @param ... Further parameters (tbd)
#' @keywords internal

faster_chans <- function(data,
                         sds = 3,
                         ...) {
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
#' @param data `eeg_epochs` object
#' @param sds standard deviations for threshold
#' @param ... Further parameters (tbd)
#' @keywords internal

faster_epochs <- function(data, sds = 3, ...) {
  chans <- channel_names(data)
  data <- data.table::as.data.table(data)
  chan_means <- data[, lapply(.SD, mean), .SDcols = chans]
  epoch_range <- data[, lapply(.SD, function(x) max(x) - min(x)),
                      .SDcols = chans,
                      by = epoch]
  epoch_range <- epoch_range[, .(Mean = rowMeans(.SD)), by = epoch]
  epoch_range <- abs(scale(epoch_range$Mean)) > sds

  epoch_diffs <- data[, lapply(.SD, mean),
                      .SDcols = chans,
                      by = epoch][, lapply(.SD, function(x) x - mean(x)),
                                  .SDcols = chans][ ,
                                                    .(Mean = rowMeans(.SD))]
  epoch_diffs <- abs(scale(epoch_diffs$Mean)) > sds

  epoch_vars <- data[, lapply(.SD, var), .SDcols = chans,
                     by = epoch][, apply(.SD, 1, mean),
                                 .SDcols = chans]
  epoch_vars <- abs(scale(epoch_vars)) > sds

  bad_epochs <- matrix(c(rowSums(epoch_vars) > 0,
                         rowSums(epoch_range) > 0,
                         rowSums(epoch_diffs) > 0),
                       ncol = 3)
  bad_epochs <- apply(bad_epochs, 1, any)
  bad_epochs
}

#' FASTER detection of bad channels in single epochs
#'
#' @param data `eeg_epochs` object.
#' @param exclude Channels to be ignored.
#' @param ... further parameters (tbd)
#' @keywords internal

faster_cine <- function(data,
                        exclude = NULL,
                        max_bad = 10,
                        ...) {

  # get xyz coords only
  xyz_coords <- data$chan_info[, c("electrode",
                                   "cart_x",
                                   "cart_y",
                                   "cart_z")]
  # check for rows with missing values
  missing_values <- apply(xyz_coords, 1,
                          function(x) any(is.na(x)))
  #remove any rows with missing values
  xyz_coords <- xyz_coords[!missing_values, ]

  keep_chans <- names(data$signals) %in% xyz_coords$electrode

  if (is.null(exclude)) {
    chan_means <- colMeans(data$signals)
  } else {
    chan_means <- colMeans(data$signals[!colnames(data$signals) %in% exclude])
  }

  epochs <- split(data$signals,
                  data$timings$epoch)
  n_epochs <- length(epochs)

  # Work out which chans are bad in each epoch according to FASTER
  bad_chans <- lapply(epochs,
                      faster_epo_stat,
                      exclude = exclude,
                      chan_means = chan_means)

  # remove channel names that are for channels we have no locations for
  bad_chans <- lapply(bad_chans,
                      function(x) x[x %in% xyz_coords$electrode])

  # remove any epochs where there were no bad channels
  bad_chans <- bad_chans[lapply(bad_chans, length) > 0]

  n_bads <- vapply(bad_chans, length, FUN.VALUE = integer(1))
  broken_epochs <- names(n_bads)[n_bads > max_bad]
  repairable_epochs <- names(n_bads)[n_bads <= max_bad]

  if (length(repairable_epochs) >= 1) {
    bad_chans <- bad_chans[repairable_epochs]
  }

  # Get a transfer matrix for each epoch
  bad_coords <- lapply(bad_chans,
                       function(x) interp_weights(xyz_coords,
                                                  x))

  bad_coords <- bad_coords[lapply(bad_coords, length) > 0]

  # If there's nothing bad in any epoch, return the data
  if (length(bad_coords) == 0) {
    return(data)
  }

  bad_epochs <- names(bad_coords)

  new_epochs <- lapply(bad_epochs,
                       function(x) interp_chans(epochs[[x]],
                                                bad_chans[[x]],
                                                !keep_chans,
                                                bad_coords[[x]]))

  epochs <- replace(epochs,
                    bad_epochs,
                    new_epochs)
  epochs <- epochs[!names(epochs) %in% broken_epochs]
  epochs <- data.table::rbindlist(epochs)
  data$signals <- tibble::as_tibble(epochs)
  data$reject$cine_list <- bad_chans
  data$reject$cine_total <- unlist(lapply(bad_chans,
                                          length))

  if (length(bad_chans) > 0) {
    message(paste0(length(bad_chans), " of ", n_epochs, " epochs had at least one channel interpolated."))
    message(paste0("Max number of channels interpolated in one epoch: ", max(data$reject$cine_total)))
  } else {
    message("No channels were interpolated in single epochs.")
  }
  data
}

#' @noRd
interp_weights <- function(xyz_coords, x) {

  xyz_coords[, c("cart_x", "cart_y", "cart_z")] <-
    norm_sphere(xyz_coords[, c("cart_x", "cart_y", "cart_z")])

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

quick_hurst <- function(data) {
  n <- nrow(data)
  data <- data.table::data.table(data)
  dat_cumsum <- data[, lapply(.SD, cumsum)]
  rs <- dat_cumsum[, lapply(.SD, max)] - dat_cumsum[, lapply(.SD, min)]
  rs <- rs / data[, lapply(.SD, stats::sd)]#column_sd
  as.numeric(log(rs) / log(n))
}

#' Calculate statistics for each channel in an epoch and identify bad channels
#'
#' @param data a matrix of signals from a single epoch
#' @keywords internal

faster_epo_stat <- function(data,
                            exclude = NULL,
                            sds = 3,
                            chan_means) {

  if (!is.null(exclude)) {
    data <- select(data, -exclude)
  }

  measures <-
    data.frame(vars = matrixStats::colVars(as.matrix(data)),
               medgrad = matrixStats::colMedians(diff(as.matrix(data))),
               range_diff = t(diff(t(matrixStats::colRanges(as.matrix(data))))),
               chan_dev = abs(colMeans(data) - chan_means)
  )
  # Check if any measure is above threshold of standard deviations
  # for some reason FASTER median centres all measures.

  bad_chans <- rowSums(abs(scale(measures)) > sds) > 0
  bad_chans <- names(data)[bad_chans]
  bad_chans
}

#' @describeIn ar_FASTER Run FASTER on `eeg_group` objects
#' @param EOG names of EOG channels to be used when computed maximum EOG values.
#' @export
ar_FASTER.eeg_group <- function(data,
                                exclude = NULL,
                                test_chans = TRUE,
                                test_epochs = TRUE,
                                test_cine = TRUE,
                                EOG = NULL,
                                ...) {

  if (!inherits(data, "eeg_evoked")) {
    stop("FASTER for grouped data only works for group ERPs.")
  }

  n_epochs <- length(unique(epochs(data)$epoch))
  n_participants <- length(unique(epochs(data)$participant_id))

  orig_names <- channel_names(data)

  data_chans <- orig_names[!(orig_names %in% data$reference$ref_chans)]
  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      exclude <- orig_names[exclude]
    }
    message("Excluding channel(s):",
            paste(exclude, ""))
    data_chans <- data_chans[!(data_chans %in% exclude)]
  }

  if (!is.null(EOG)) {
    if (!all(EOG %in% orig_names)) {
      stop("EOG channels not specified correctly, please double-check.")
    }
  }

  all_data <- as.data.frame(data)

  participant_mean <-
    dplyr::ungroup(all_data) %>%
    dplyr::group_by(participant_id) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(data_chans),
        mean)
      )

  channel_mean <-
    colMeans(all_data[, data_chans])

  part_diffs <-
    abs(sweep(participant_mean[, data_chans],
              2,
              channel_mean, "-"))

  participant_var <-
    dplyr::ungroup(all_data) %>%
    dplyr::group_by(participant_id) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(data_chans),
        var)
    )

  channel_vars <-
    matrixStats::colVars(as.matrix(all_data[, data_chans]))

  part_var_diffs <-
    abs(sweep(participant_var[, data_chans],
              2,
              channel_vars, "-"))

  part_ranges <-
    dplyr::ungroup(all_data) %>%
    dplyr::group_by(participant_id) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(data_chans),
        ~diff(range(.)))
    )

  scaled_diffs <- abs(scale(part_diffs[, data_chans]))
  scaled_vars <- abs(scale(part_var_diffs[, data_chans]))
  scaled_ranges <- abs(scale(part_ranges[, data_chans]))

  # to be added - checks of residual HEOG/VEOG

  if (!is.null(EOG)) {
    participant_eog <-
      dplyr::ungroup(all_data) %>%
      dplyr::group_by(participant_id) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(EOG),
          max)
      )
    scaled_eog <- abs(scale(matrixStats::rowMaxs(as.matrix(participant_eog[, EOG]))))
    data.frame(
      participant_id = unique(all_data$participant_id),
      deviance = rowSums(scaled_diffs > 3),
      variance = rowSums(scaled_vars > 3),
      range = rowSums(scaled_ranges > 3),
      eog = rowSums(scaled_eog > 3)
    )
  } else {
    data.frame(
      participant_id = unique(all_data$participant_id),
      deviance = rowSums(scaled_diffs > 3),
      variance = rowSums(scaled_vars > 3),
      range = rowSums(scaled_ranges > 3)
      )
  }

}
