#' Baseline correction
#'
#' Used to remove the mean of a specified time period from the data. Currently
#' only performs subtractive baseline. With a data frame, searches for
#' "electrode" and "epoch" columns, and groups on these when found. An electrode
#' column is always required; an epoch column is not.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to be baseline corrected.
#' @param time_lim Numeric character vector (e.g. time_lim <- c(-.1, 0)). If
#'   none given, defaults to mean of the whole of each epoch if the data is
#'   epoched, or the channel mean if the data is continuous.
#' @param verbose Defaults to TRUE. Output descriptive messages to console.
#' @param ... other parameters to be passed to functions
#' @examples
#' rm_baseline(demo_epochs)
#' rm_baseline(demo_epochs, c(-.1, 0))
#' @export

rm_baseline <- function(data,
                        time_lim = NULL,
                        ...) {
  UseMethod("rm_baseline", data)
}


#' @describeIn rm_baseline remove baseline from continuous \code{eeg_data}
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
  } else {
    base_times <- select_times(data,
                               time_lim = time_lim)
    baseline_dat <- colMeans(base_times$signals)
    if (verbose) {
      message(paste("Baseline:", time_lim, "s"))
    }
    data$signals <- as.matrix(data$signals)
  }

  data$signals <- baseline_cont(data$signals,
                                baseline_dat)
  data$signals <- tibble::as_tibble(data$signals)
  data
}

#' @describeIn rm_baseline Remove baseline from eeg_epochs
#' @export

rm_baseline.eeg_epochs <- function(data,
                                   time_lim = NULL,
                                   verbose = TRUE,
                                   ...) {

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

  } else {
    base_times <- select_times(data,
                               time_lim = time_lim)
    base_times$signals <- as.matrix(base_times$signals)
    n_bl_times <- length(unique(base_times$timings$time))
    dim(base_times$signals) <- c(n_bl_times, n_epochs, n_chans)
    base_times <- colMeans(base_times$signals)

    data$signals <- as.matrix(data$signals)
    dim(data$signals) <- c(n_times, n_epochs, n_chans)
    data$signals <- baseline_epo(data$signals, base_times)
  }
  #Reshape and turn back into data frame
  data$signals <- array(data$signals,
                       dim = c(n_epochs * n_times, n_chans))
  data$signals <- tibble::as_tibble(data$signals)

  names(data$signals) <- elecs
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
  } else{
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
#' @describeIn rm_baseline Method for \code{eeg_tfr} objects
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

  if (length(dim(data$signals)) == 4) {
    epoched <- TRUE
  } else {
    epoched <- FALSE
  }

  bline <- select_times(data, time_lim)
  if (epoched) {
    #bline <- colMeans(colMeans(bline$signals, na.rm = TRUE), na.rm = TRUE)
    bline <- apply(bline$signals, c(1, 3, 4), mean, na.rm = TRUE)
  } else {
    bline <- colMeans(bline$signals, na.rm = TRUE)
  }
  # This function implements the various baseline correction types
  do_corrs <- function(data,
                       type,
                       bline) {
    switch(type,
           "divide" = ((data - bline) / bline) * 100,
           "pc" = ((data - bline) / bline) * 100 - 100,
           "absolute" = data - bline,
           "db" = 10 * log10(data / bline),
           "ratio" = data / bline
    )
  }

  orig_dims <- dim(data$signals)

  orig_dimnames <- dimnames(data$signals)

  if (epoched) {

    data$signals <- sweep(data$signals,
                          c(1, 3, 4),
                          bline,
                          do_corrs,
                          type = type)


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

#' @describeIn rm_baseline Method for \code{eeg_evoked} objects
#' @export
rm_baseline.eeg_evoked <- function(data,
                                   time_lim = NULL,
                                   verbose = TRUE,
                                   ...) {

  # Edit to handle cases where multiple epochs/participants, now that the internal structure differs
  if (is.null(time_lim)) {
    data$signals <- data.table::as.data.table(as.data.frame(data))

  } else {
    orig_cols <- channel_names(data)
    data$signals <- data.table::as.data.table(as.data.frame(data))
    data$signals <-
      data$signals[, c(orig_cols) := lapply(.SD,
                                            function(x) x - mean(x[time > time_lim[1] & time < time_lim[2]])),
       .SDcols = orig_cols,
       by = epoch]
    data$signals <- tibble::as_tibble(data$signals[, ..orig_cols])

  }
  data
}
