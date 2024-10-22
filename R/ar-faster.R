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
#' @references Nolan, Whelan & Reilly (2010). FASTER: Fully Automated
#' Statistical Thresholding for EEG artifact Rejection. J Neurosci Methods.
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

  data <- validate_input(data)

  metadata <- prepare_data(data = data, exclude = exclude)

  if (test_chans) {
    data <- handle_bad_channels(data = data, data_chans = metadata$data_chans)
  } else {
    message("Skipping global bad channel detection...")
  }

  if (test_epochs) {
    data <- handle_bad_epochs(data)
  } else {
    message("Skipping globally bad epoch detection...")
  }

  if (test_cine) {
    data <- handle_bad_channels_in_epochs(data = data, exclude = exclude)
  } else {
    message("Skipping detection of locally bad channels in epochs...")
  }

  data <- restore_original_reference(data = data, metadata = metadata)

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
  chan_vars <- scale(apply(data, 2, stats::var))
  chan_corrs <- scale(colMeans(abs(stats::cor(data))))

  # Optimized version
  bad_chans <- (abs(chan_hurst) > sds) |
               (abs(chan_vars) > sds) |
               (abs(chan_corrs) > sds)

  bad_chans
}

#' Perform global bad epoch detection for FASTER
#'
#' @param data `eeg_epochs` object
#' @param sds standard deviations for threshold
#' @param ... Further parameters (tbd)
#' @keywords internal

faster_epochs <- function(data, sds = 3, ...) {
  chan_means <- colMeans(data$signals)
  data$signals <- split(data$signals, data$timings$epoch)
  data$signals <- lapply(data$signals, as.matrix)
  epoch_ranges <- lapply(data$signals,
                         function(x) {
                           matrixStats::rowDiffs(
                             matrixStats::colRanges(x)
                           )
                         })
  epoch_ranges <- matrix(unlist(epoch_ranges),
                         ncol = length(epoch_ranges))
  epoch_ranges <- colMeans(epoch_ranges)

  epoch_diffs <- lapply(data$signals,
                        colMeans)
  epoch_diffs <- matrix(unlist(epoch_diffs),
                        ncol = length(epoch_diffs))
  epoch_diffs <- epoch_diffs - chan_means
  epoch_diffs <- colMeans(abs(epoch_diffs))

  epoch_vars <- lapply(data$signals,
                       function(x) matrixStats::colVars(x))
  epoch_vars <- matrix(unlist(epoch_vars),
                       ncol = length(epoch_vars))
  epoch_vars <- colMeans(epoch_vars)

  measures <- matrix(c(epoch_ranges,
                       epoch_diffs,
                       epoch_vars),
                     ncol = 3)

  measures <- abs(scale(measures)) >= sds
  measures <- rowSums(measures) > 0
  measures
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

  # Get xyz coords only
  xyz_coords <- data$chan_info[, c("electrode", "cart_x", "cart_y", "cart_z")]

  # Check for rows with missing values (vectorized operation)
  missing_values <- rowSums(is.na(xyz_coords)) > 0

  # Remove any rows with missing values
  xyz_coords <- xyz_coords[!missing_values, ]

  # Use %in% for faster matching
  keep_chans <- names(data$signals) %in% xyz_coords$electrode

  # Precompute chan_means
  if (is.null(exclude)) {
    chan_means <- colMeans(data$signals)
  } else {
    chan_means <- colMeans(data$signals[, !colnames(data$signals) %in% exclude, drop = FALSE])
  }

  # Use data.table for faster operations
  epochs <- data.table::as.data.table(data$signals)
  epochs[, epoch := data$timings$epoch]

  n_epochs <- uniqueN(epochs$epoch)

  # Work out which chans are bad in each epoch according to FASTER
  bad_chans <- epochs[, .(bad_chan = faster_epo_stat(.SD, exclude, chan_means = chan_means)), by = epoch]

  # Filter bad channels
  bad_chans <- bad_chans[bad_chan %in% xyz_coords$electrode]

  # Count bad channels per epoch
  n_bads <- bad_chans[, .N, by = epoch]

  broken_epochs <- n_bads[N > max_bad, epoch]
  repairable_epochs <- n_bads[N <= max_bad & N > 0, epoch]

  if (length(repairable_epochs) >= 1) {
    bad_chans <- bad_chans[epoch %in% repairable_epochs]
  }

  # Get a transfer matrix for each epoch
  bad_coords <- lapply(bad_chans,
                       function(x) interp_weights(xyz_coords,
                                                  x))

  bad_coords <- bad_coords[lapply(bad_coords, length) > 0]

  # If there's nothing bad in any epoch, return the data
  if (nrow(bad_coords) == 0) {
    return(data)
  }

  bad_epochs <- bad_coords$epoch

  # Use data.table for faster operations
  new_epochs <- epochs[epoch %in% bad_epochs,
                       interp_chans(.SD,
                                    bad_chans[epoch == .BY$epoch, bad_chan],
                                    !keep_chans,
                                    bad_coords[epoch == .BY$epoch, coords[[1]]]),
                       by = epoch]

  # Update epochs
  epochs <- epochs[!epoch %in% broken_epochs]
  epochs[epoch %in% bad_epochs, names(epochs) := new_epochs]

  data$signals <- as.data.frame(epochs[, !c("epoch")])
  data$reject$cine_list <- bad_chans[, .(bad_chan = list(bad_chan)), by = epoch]
  data$reject$cine_total <- bad_chans[, .N, by = epoch]$N

  if (nrow(bad_chans) > 0) {
    message(sprintf("%d of %d epochs had at least one channel interpolated.",
                    uniqueN(bad_chans$epoch), n_epochs))
    message(sprintf("Max number of channels interpolated in one epoch: %d",
                    max(data$reject$cine_total)))
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
        mean,
        na.rm = TRUE)
    )

  channel_mean <-
    colMeans(all_data[, data_chans], na.rm = TRUE)

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
        var,
        na.rm = TRUE)
    )

  channel_vars <-
    matrixStats::colVars(as.matrix(all_data[, data_chans]), na.rm = TRUE)

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
        ~diff(range(., na.rm = TRUE)))
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
          max, na.rm = TRUE)
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

#' Validate input data for FASTER
#'
#' This function checks if the input data is valid and has the required channel information.
#'
#' @param data An object of class `eeg_epochs`
#' @return The validated data object
#' @keywords internal
validate_input <- function(data) {
  check_ci_str(data$chan_info)

  if (is.null(channels(data))) {
    stop("No channel information found, ar_FASTER() requires (at least some) channel locations.")
  } else {
    channels(data) <- validate_channels(channels(data),
                                        channel_names(data))
  }
  data
}

#' Prepare data for FASTER processing
#'
#' This function prepares the data for processing, handling exclusions and storing original reference information.
#'
#' @param data An object of class `eeg_epochs`
#' @param exclude Channels to be ignored by FASTER
#' @return A list containing prepared data information
#' @keywords internal
prepare_data <- function(data, exclude) {
  orig_ref <- data$reference$ref_chans
  excluded <- data$reference$excluded

  orig_chan_info <- channels(data)
  orig_names <- channel_names(data)
  data_chans <- orig_names[!(orig_names %in% data$reference$ref_chans)]

  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      exclude <- orig_names[exclude]
    }
    message("Excluding channel(s):", paste(exclude, ""))
    data_chans <- data_chans[!(data_chans %in% exclude)]
  }

  list(orig_ref = orig_ref,
       orig_chan_info = orig_chan_info,
       data_chans = data_chans)
}

#' Handle bad channels globally
#'
#' This function detects and handles globally bad channels using the FASTER algorithm.
#'
#' @param data An object of class `eeg_epochs`
#' @param data_chans A vector of channel names to be processed
#' @return The data object with bad channels handled
#' @keywords internal
handle_bad_channels <- function(data, data_chans) {
  bad_chans <- faster_chans(data$signals[, data_chans])
  bad_chan_n <- data_chans[bad_chans]
  message(paste("Globally bad channels:",
                paste(bad_chan_n, collapse = " ")))

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
  data
}

#' Handle bad epochs globally
#'
#' This function detects and handles globally bad epochs using the FASTER algorithm.
#'
#' @param data An object of class `eeg_epochs`
#' @return The data object with bad epochs handled
#' @keywords internal
handle_bad_epochs <- function(data) {
  bad_epochs <- faster_epochs(data)
  bad_epochs <- unique(data$timings$epoch)[bad_epochs]
  message(paste("Globally bad epochs:",
                paste(bad_epochs,
                      collapse = " ")))
  data$reject$bad_epochs <- bad_epochs
  data <- select_epochs(data,
                        epoch_no = bad_epochs,
                        keep = FALSE)
  data
}

#' Handle bad channels in individual epochs
#'
#' This function detects and handles bad channels in individual epochs using the FASTER algorithm.
#'
#' @param data An object of class `eeg_epochs`
#' @param exclude Channels to be ignored by FASTER
#' @return The data object with bad channels in individual epochs handled
#' @keywords internal
handle_bad_channels_in_epochs <- function(data, exclude) {
  data <- faster_cine(data, exclude)
  data
}

#' Restore original reference
#'
#' This function restores the original reference if it existed.
#'
#' @param data An object of class `eeg_epochs`
#' @param metadata A list containing metadata information, including the original reference
#' @return The data object with the original reference restored
#' @keywords internal
restore_original_reference <- function(data, metadata) {
  if (!is.null(metadata$orig_ref)) {
    data <- eeg_reference(data,
                          ref_chans = metadata$orig_ref,
                          exclude = metadata$excluded,
                          verbose = FALSE)
  }
  data
}

