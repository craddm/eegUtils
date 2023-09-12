#' Baseline correction
#'
#' Used to correct data using the mean of a specified time period. For
#' time-domain data, this will subtract the mean from all data. For `eeg_tfr`
#' objects, a variety of methods are available, including subtraction, and
#' conversion to "dB" change. With a data frame, it will search for "electrode"
#' and "epoch" columns, and groups on these when found. An electrode column is
#' always required; an epoch column is not. Note that baseline correction is
#' always applied on single-trial basis. For baseline correction based on
#' subtraction, this makes no difference compared to averaging first and then
#' baseline correcting, but for divisive measures used with time-frequency data,
#' this distinction can be very important, and can lead to counterintuitive
#' results.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to be baseline corrected.
#' @param time_lim Numeric character vector (e.g. time_lim <- c(-.1, 0))
#'   defining the time period to use as a baseline. If the value is NULL, it
#'   uses the mean of the whole of each epoch if the data is epoched, or the
#'   channel mean if the data is continuous.
#' @param verbose Defaults to TRUE. Output descriptive messages to console.
#' @param ... other parameters to be passed to functions
#' @return An `eegUtils` object or a `data.frame`, depending on the input.
#' @examples
#' rm_baseline(demo_epochs)
#' rm_baseline(demo_epochs, c(-.1, 0))
#' @export

rm_baseline <- function(data,
                        time_lim = NULL,
                        ...) {
  UseMethod("rm_baseline", data)
}


#' @describeIn rm_baseline remove baseline from continuous `eeg_data`
#' @export

rm_baseline.eeg_data <- function(data,
                                 time_lim = NULL,
                                 verbose = TRUE,
                                 ...) {

  if (is.null(time_lim)) {
    data$signals <- as.matrix(data$signals)
    baseline_dat <- colMeans(data$signals)
    if (verbose) {
      message("Removing channel means...")
    }
    data$baseline <- range(data$timings$time)
  } else {
    base_times <- select_times(data,
                               time_lim = time_lim)
    baseline_dat <- colMeans(base_times$signals)
    if (verbose) {
      message(paste0(
        "Baseline: ",
        time_lim[1],
        " - ",
        time_lim[2],
        "s")
        )
    }
    data$signals <- as.matrix(data$signals)
    data$baseline <- time_lim
  }

  data$signals <- baseline_cont(data$signals,
                                baseline_dat)
  data$signals <- tibble::as_tibble(data$signals)

  data
}

#' @describeIn rm_baseline Remove baseline from `eeg_epochs`
#' @export

rm_baseline.eeg_epochs <- function(data,
                                   time_lim = NULL,
                                   verbose = TRUE,
                                   ...) {

  # modify to handle group objects
  n_epochs <- length(unique(data$timings$epoch))
  n_times <- length(unique(data$timings$time))
  n_chans <- ncol(data$signals)
  elecs <- names(data$signals)

  # I calculate the baseline for each epoch and subtract it; this makes no
  # difference to ERP later, but centres each epoch on zero.
  if (is.null(time_lim)) {

    if (verbose) {
      message("Removing channel means per epoch...")
    }
    # reshape to 3D matrix
    orig_chans <- channel_names(data)

    data$signals <- as.matrix(data$signals)
    dim(data$signals) <- c(n_times, n_epochs, n_chans)
    # colMeans gives an n_epochs * n_channels matrix - i.e. baseline value for
    # each epoch and channel
    baseline_dat <- colMeans(data$signals)
    # now we go through each timepoint subtracting the baseline values
    data$signals <- baseline_epo(data$signals,
                                 baseline_dat)
    data$baseline <- range(data$timings$time)
  } else {
    if (verbose) {
      message(paste0(
        "Baseline: ",
        time_lim[1],
        " - ",
        time_lim[2],
        "s"))
    }
    base_times <- get_epoch_baselines(data,
                                      time_lim)

    data$signals <- as.matrix(data$signals)
    dim(data$signals) <- c(n_times, n_epochs, n_chans)
    data$signals <- baseline_epo(data$signals, base_times)
    data$baseline <- time_lim
  }
  #Reshape and turn back into data frame
  data$signals <- array(data$signals,
                        dim = c(n_epochs * n_times,
                                n_chans))
  colnames(data$signals) <- elecs
  data$signals <- tibble::as_tibble(data$signals)
  data
}

#' @describeIn rm_baseline Legacy method for data.frames
#' @export
rm_baseline.data.frame <- function(data,
                                   time_lim = NULL,
                                   verbose = TRUE,
                                   ...) {

  warning("rm_baseline.data.frame will be deprecated.")

  if (!("time" %in% colnames(data))) {
    stop("Time dimension is required.")
  }

  if (length(time_lim) == 1) {
    stop("time_lim should specify the full time range.")
  }

  # if the data is epoched, group by electrode and epoch; otherwise, just by
  # electrode.

  if ("epoch" %in% colnames(data)) {
    data <- dplyr::group_by(data,
                            electrode,
                            epoch,
                            add = TRUE)
  } else {
    data <- dplyr::group_by(data,
                            electrode,
                            add = TRUE)
  }

  if (is.null(time_lim)) {
    # if no time_lim provided, just delete mean of all time points
    data <- dplyr::mutate(data,
                          amplitude = amplitude - mean(amplitude))
  } else {

    data_sel <- dplyr::filter(data,
                              time >= time_lim[1],
                              time <= time_lim[2])
    baseline <- dplyr::summarise(data_sel,
                                 bl = mean(amplitude))
    # This is relatively memory intensive - not so bad now but would prefer
    # another way. Could get extremely painful with time-frequency data.
    data <- dplyr::left_join(data,
                             baseline)
    data <- dplyr::mutate(data,
                          amplitude = amplitude - bl)
    data <- dplyr::select(data,
                          -bl)
  }
  data <- ungroup(data)
  data
}

#' @param type Type of baseline correction to apply. Options are ("divide",
#'   "ratio", "absolute", "db", and "pc")
#' @describeIn rm_baseline Method for `eeg_tfr` objects
#' @export
rm_baseline.eeg_tfr <- function(data,
                                time_lim = NULL,
                                type = "divide",
                                verbose = TRUE,
                                ...) {

  valid_types <- c("absolute",
                   "divide",
                   "pc",
                   "ratio",
                   "db")

  if (!(type %in% valid_types)) {
    stop("Unknown baseline type ", type)
  }

  is_group_tfr <- inherits(data,
                           "eeg_group")

  orig_dims <- dimnames(data$signals)

  if ("epoch" %in% names(orig_dims)) {
    epoched <- TRUE
  } else {
    epoched <- FALSE
  }

  if (!is.null(time_lim)) {
    if (verbose) {
      message(paste0("Baseline: ",
                     time_lim[1],
                     "-",
                     time_lim[2],
                     " s"))
    }
    bline <- select_times(data,
                          time_lim)
    no_bline <- FALSE
  } else {
    no_bline <- TRUE
    if (verbose) {
      message(paste("Using whole epoch as baseline."))
    }
  }

  if (epoched) {

    if (is_group_tfr) {
      if (no_bline) {
        bline <- apply(data$signals,
                       c(1, 3, 4, 5),
                       mean,
                       na.rm = TRUE)
      } else {
        bline <- apply(bline$signals,
                       c(1, 3, 4, 5),
                       mean,
                       na.rm = TRUE)
      }

    } else {
      if (no_bline) {
        bline <- apply(data$signals,
                       c(1, 3, 4),
                       mean,
                       na.rm = TRUE)
      } else {
        bline <- apply(bline$signals,
                       c(1, 3, 4),
                       mean,
                       na.rm = TRUE)
      }
    }
  } else {
    bline <- colMeans(bline$signals,
                      na.rm = TRUE)
  }

  # This function implements the various baseline correction types
  do_corrs <- function(data,
                       type,
                       bline) {
    switch(
      type,
      "divide" = data / bline,
      "pc" = ((data - bline) / bline) * 100,
      "absolute" = data - bline,
      "db" = 10 * log10(data / bline),
      "ratio" = data / bline
    )
  }

  orig_dims <- dim(data$signals)

  orig_dimnames <- dimnames(data$signals)

  if (epoched) {
    if (is_group_tfr) {
      data$signals <- sweep(data$signals,
                            c(1, 3, 4, 5),
                            bline,
                            do_corrs,
                            type = type)
    } else {
      data$signals <- sweep(data$signals,
                            c(1, 3, 4),
                            bline,
                            do_corrs,
                            type = type)
    }
  } else {
    data$signals <- apply(data$signals,
                          1,
                          do_corrs,
                          type = type,
                          bline = bline)
    data$signals <- aperm(data$signals,
                          c(2, 1))
  }

  dim(data$signals) <- orig_dims

  dimnames(data$signals) <- orig_dimnames
  data$freq_info$baseline <- type
  data$freq_info$baseline_time <- time_lim
  data
}

#' @describeIn rm_baseline Method for `eeg_evoked` objects
#' @export
rm_baseline.eeg_evoked <- function(data,
                                   time_lim = NULL,
                                   verbose = TRUE,
                                   ...) {

  orig_cols <- channel_names(data)
  n_times <- length(unique(data$timings$time))
  n_epochs <- nrow(unique(epochs(data)[, c("epoch", "participant_id")]))
  n_participants <- length(unique(epochs(data)$participant_id))
  n_chans <- length(orig_cols)

  if (is.null(time_lim)) {
    if (verbose) {
      message("Removing channel means...")
    }
    data$baseline <- range(data$timings$time)
  } else {
    if (verbose) {
      message(paste0(
        "Baseline: ",
        time_lim[1],
        " - ",
        time_lim[2],
        "s"))
    }
    data$baseline <- time_lim
  }

  base_times <- get_epoch_baselines(data,
                                    time_lim)

  data$signals <- as.matrix(data$signals)

  dim(data$signals) <- c(n_times,
                         n_epochs,
                         n_chans)

  data$signals <- baseline_epo(data$signals, base_times)

  data$signals <- array(data$signals,
                        dim = c(n_epochs * n_times, n_chans))
  colnames(data$signals) <- orig_cols
  data$signals <- tibble::as_tibble(data$signals)
  data
}

#' @export
rm_baseline.eeg_group <- function(data, ...) {
  message("Baseline correction support for `eeg_group` objects is currently experimental, use at own risk...!")
  NextMethod("rm_baseline")

}


#' Get epoch baselines
#'
#' Gets the baseline values for every epoch separately
#'
#' @param data data for which to calculate the baselines
#' @param time_lim time limits of the baseline period. numeric vector of length
#'   two, c(start, end)
#' @return A numeric matrix of n_epochs x n_channels.
#' @keywords internal
get_epoch_baselines <- function(data,
                                time_lim) {

  #n_epochs <- nrow(epochs(data))
  n_epochs <- nrow(unique(epochs(data)[, c("epoch", "participant_id")]))
  n_participants <- length(unique(epochs(data)$participant_id))

  n_chans <- length(channel_names(data))
  chan_names <- colnames(data$signals)

  if (is.null(time_lim)) {
    data$signals <- as.matrix(data$signals)
    n_times <- length(unique(data$timings$time))
    dim(data$signals) <- c(n_times,
                           n_epochs,
                           n_chans)
    base_times <- colMeans(data$signals)
  } else {
    base_times <- select_times(data,
                               time_lim = time_lim)
    base_times$signals <- as.matrix(base_times$signals)
    n_bl_times <- length(unique(base_times$timings$time))
    dim(base_times$signals) <- c(n_bl_times,
                                   n_epochs,
                                   n_chans)
    base_times <- colMeans(base_times$signals)
  }
  colnames(base_times) <- chan_names
  base_times
}
