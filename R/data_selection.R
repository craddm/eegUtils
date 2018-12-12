#' Select timerange
#'
#' Generic function for selecting specific time ranges from a given dataset.
#' Input can be a dataframe, or an object of class \code{eeg_data},
#' \code{eeg_epochs}, or \code{eeg_evoked}. Note this finds the closest times to
#' those specified, so the time range returned may be slightly longer or shorter
#' than that requested.
#'
#' @examples
#' ## Select timepoints from -.1 to .3
#' demo_epochs
#' short_epochs <- select_times(demo_epochs, time_lim = c(-.1, .3))
#' short_epochs
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @param data Data from which to select
#' @param ... Further arguments passed to or from other methods.
#' @family Data selection functions
#' @seealso \code{\link{select_elecs}} and \code{\link{select_epochs}}
#' @export

select_times <- function(data, ...) {
  UseMethod("select_times", data)
}

#' @param time_lim A character vector of two numbers indicating the time range
#'   to be selected e.g. c(min, max)
#' @return Data frame with only data from within the specified range.
#' @export
#' @describeIn select_times Default select times function

select_times.default <- function(data,
                                 time_lim = NULL,
                                 ...) {

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      stop("Must enter two timepoints when selecting a time range.")
    } else if (length(time_lim) == 2) {
      data <- data[data$time > time_lim[1] &
                     data$time < time_lim[2], ]
    }
  } else {
    warning("No time column found.")
  }
  data
}

#' @param df_out Returns a data frame rather than an object of the same type
#'   that was passed in.
#' @export
#' @return \code{eeg_data} object
#' @describeIn select_times Select times from an eeg_data object

select_times.eeg_data <- function(data,
                                  time_lim = NULL,
                                  df_out = FALSE,
                                  ...) {

  #data$signals <- as.data.frame(data)
  keep_rows <- find_times(data$timings, time_lim)
  data$signals <- data$signals[keep_rows, ]
  data$timings <- data$timings[keep_rows, ]
  event_rows <- data$events$event_time > time_lim[1] &
        data$events$event_time < time_lim[2]
  data$events <- data$events[event_rows, ]

  if (df_out) {
    return(as.data.frame(data))
  }

  data
}

#' @export
#' @describeIn select_times Select times in \code{eeg_epochs} objects
select_times.eeg_epochs <- function(data,
                                    time_lim,
                                    df_out = FALSE,
                                    ...) {

  keep_rows <- find_times(data$timings,
                         time_lim)

  data$signals <- data$signals[keep_rows, ]
  data$timings <- data$timings[keep_rows, ]
  event_rows <- data$events$time > time_lim[1] &
    data$events$time < time_lim[2]
  data$events <- data$events[event_rows, ]
  if (df_out) {
    return(as.data.frame(data))
  }
  data
}

#' @export
#' @describeIn select_times Select times in \code{eeg_evoked} objects
select_times.eeg_evoked <- function(data,
                                    time_lim,
                                    df_out = FALSE,
                                    ...) {

  #data$signals <- as.data.frame(data)

  keep_rows <- find_times(data$timings,
                          time_lim)

  data$signals <- data$signals[keep_rows, ]
  data$timings <- data$timings[keep_rows, ]

  if (df_out) {
    return(data$signals)
  }
  #data$signals$time <- NULL
  data
}

#' @describeIn select_times Select times from an \code{eeg_tfr} object
#' @export
select_times.eeg_tfr <- function(data,
                                 time_lim = NULL,
                                 df_out = FALSE,
                                 ...){

  keep_rows <- find_times(data$timings, time_lim)
  data$timings <- data$timings[keep_rows, ]
  if (length(data$dimensions) == 3) {
    data$signals <- data$signals[keep_rows, , , drop = FALSE]
  } else {
    data$signals <- data$signals[keep_rows, , , , drop = FALSE]
  }
  data
}

#' Find times in an eeg_* object
#'
#' Internal function to find the rows corresponding to the selected time limits
#'
#' @param timings timing information from the EEG object.
#' @param time_lim character vector of the time limits
#' @return logical vector of selected timepoints
#' @keywords internal
find_times <- function(timings,
                       time_lim) {

  if (length(time_lim) == 2) {
    keep_rows <- timings$time > time_lim[1] & timings$time < time_lim[2]
  } else {
  warning("Must enter two timepoints when selecting a time range;
          using whole range.")
    keep_rows <- rep(TRUE, length = length(timings$time))
  }
  keep_rows
}

#' Select electrodes from a given dataset
#'
#' This is a generic function for selection of electrodes from an EEG dataset.
#'
#' @examples
#' names(demo_epochs$signals)
#' keep_A5 <- select_elecs(demo_epochs, electrode = "A5")
#' remove_A5 <- select_elecs(demo_epochs, electrode = "A5", keep = FALSE)
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @param data An EEG dataset.
#' @param ... Arguments used with related methods
#' @family Data selection functions
#' @seealso \code{\link{select_times}} and \code{\link{select_epochs}}
#' @export

select_elecs <- function(data, ...) {
  UseMethod("select_elecs", data)
}

#' @param electrode A character vector of electrode labels for selection or
#'   removal.
#' @param keep Defaults to TRUE. Set to false to *remove* the selected
#'   electrodes.
#' @return Data frame with only data from the chosen electrodes
#' @describeIn select_elecs Select electrodes from a generic data frame.
#' @export

select_elecs.default <- function(data,
                                 electrode = NULL,
                                 keep = TRUE,
                                 ...) {

  if ("electrode" %in% names(data)) {
    if (all(electrode %in% data$electrode)) {
      if (keep) {
        data <- data[data$electrode %in% electrode, ]
      } else {
        data <- data[!data$electrode %in% electrode, ]
      }
    } else {
      warning(paste("Electrode(s) not found:",
                    electrode[!electrode %in% data$electrode],
                    ". Returning all data."))
    }
  } else {
    if (all(electrode %in% names(data))) {
      if (keep) {
        data <- data[, names(data) %in% electrode, drop = FALSE]
      } else {
        data <- data[, !names(data) %in% electrode, drop = FALSE]
      }
    }
  }
  data
}

#' @param df_out Defaults to FALSE. Set to TRUE to return a dataframe rather
#'   than an \code{eeg_data} object.
#' @return \code{eeg_data} object with selected electrodes removed/kept.
#' @export
#' @describeIn select_elecs Select electrodes from a \code{eeg_data} object.

select_elecs.eeg_data <- function(data,
                                  electrode,
                                  keep = TRUE,
                                  df_out = FALSE, ...) {

  if (all(electrode %in% names(data$signals))) {

    if (keep) {
      data$signals <- data$signals[colnames(data$signals) %in% electrode]
    } else {
      data$signals <- data$signals[!colnames(data$signals) %in% electrode]
    }

    if (!is.null(data$chan_info)) {
      data$chan_info <- data$chan_info[data$chan_info$electrode %in% names(data$signals), ]
    }

  } else {
    warning("Electrode(s) not found:",
        electrode[!electrode %in% colnames(data$signals)],
        ". Returning all data.")
    return(data)
  }

  if (df_out) {
    return(as.data.frame(data))
  }
  data
}

#' @describeIn select_elecs Select electrode from an eeg_evoked object
#' @export
select_elecs.eeg_evoked <- function(data,
                                    electrode = NULL,
                                    keep = TRUE,
                                    df_out = FALSE,
                                    ...) {

  sig_names <- electrode %in% names(data$signals)

  if (!all(sig_names)) {
    warning("Electrode(s) not found:",
            electrode[!electrode %in% names(data$signals)],
            ". Returning all data.")
    return(data)
  }

  sig_names <- names(data$signals) %in% electrode

  if (!keep) {
    sig_names <- !sig_names
  }

  data$signals <- data$signals[, sig_names, drop = FALSE]

  if (!is.null(data$chan_info)) {
    data$chan_info <- data$chan_info[data$chan_info$electrode %in% names(data$signals), ]
  }

  if (df_out) {
    return(as.data.frame(data))
  }
  data
}

#' @param component Component to select
#' @describeIn select_elecs Select components from \code{eeg_ICA} objects.
#' @export
select_elecs.eeg_ICA <- function(data,
                                 component,
                                 keep = TRUE,
                                 df_out = FALSE,
                                 ...) {

  if (!all(component %in% names(data$signals))) {
    stop("Component(s) ", component, " not found.")
  }

  comps <- names(data$signals) %in% component
  if (!keep) {
    comps <- !comps
  }
  data$mixing_matrix <- data$mixing_matrix[,
                                           c(comps, TRUE),
                                           drop = FALSE]
  data$unmixing_matrix <- data$unmixing_matrix[,
                                             c(comps,
                                               TRUE),
                                             drop = FALSE]
  data$signals <- data$signals[,
                               comps,
                               drop = FALSE]
  data
}

#'@importFrom abind asub
#'@describeIn select_elecs Select electrodes from \code{eeg_tfr} objects.
#'@export
select_elecs.eeg_tfr <- function(data,
                                 electrode,
                                 keep = TRUE,
                                 df_out = FALSE,
                                 ...) {

  elec_dim <- which(data$dimensions == "electrode")
  data_elecs <- dimnames(data$signals)[[elec_dim]]
  sig_names <- electrode %in% data_elecs

  if (!all(sig_names)) {
    warning("Electrode(s) not found:",
            electrode[!electrode %in% names(data$signals)],
            ". Returning all data.")
    return(data)
  }

  sig_names <- data_elecs %in% electrode

  if (!keep) {
    sig_names <- !sig_names
  }

  data$signals <- abind::asub(data$signals,
                              sig_names,
                              elec_dim,
                              drop = FALSE)

  if (!is.null(data$chan_info)) {
    data$chan_info <- data$chan_info[data$chan_info$electrode %in% data_elecs[sig_names], ]
  }

  if (df_out) {
    return(as.data.frame(data))
  }

  data
}

#' Select epochs
#'
#' This is a generic function for selecting epochs from an epoched data set.
#'
#' @examples
#' select_epochs(demo_epochs, epoch_no = 1:5)
#' (demo_ica <- run_ICA(demo_epochs))
#' select_epochs(demo_ica, epoch_no = 1:5)
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @param data \code{eeg_epochs} object from which to select epochs.
#' @param ... Parameters passed to specific methods
#' @family data selection functions
#' @seealso \code{\link{select_times}} and \code{\link{select_elecs}}
#' @export

select_epochs <- function(data, ...) {
  UseMethod("select_epochs", data)
}

#' @describeIn select_epochs Select from generic object
#' @export

select_epochs.default <- function(data, ...) {

  warning(paste("select_epochs does not know how to handle object of class",
                class(data),
                "and can only be used on eeg_epochs objects."))
}

#' @describeIn select_epochs Select epochs from \code{eeg_data} objects.
#' @export

select_epochs.eeg_data <- function(data, ...) {
  if (data$continuous) {
    stop("Data is not epoched.")
  } else {
    warning("oops, shouldn't end up here.")
  }
}

#' @param epoch_events Select epochs containing any of the specified events. Can
#'   be numeric or a character string. Will override any epoch_no input.
#' @param epoch_no Select epochs by epoch number.
#' @param keep Defaults to TRUE, meaning select the specified epochs. Set to
#'   FALSE to remove specified epochs.
#' @param df_out Output a data.frame instead of an eeg_data object.
#' @describeIn select_epochs Selection of epochs from \code{eeg_epochs} objects.
#' @export

select_epochs.eeg_epochs <- function(data,
                                     epoch_events = NULL,
                                     epoch_no = NULL,
                                     keep = TRUE,
                                     df_out = FALSE,
                                     ...) {

  # First check if epoch_events has been passed; if it's numeric, select epochs
  # based on event_type. If it's a character vector, check if those labels exist
  # in the data.

  if (!is.null(epoch_events)) {
    epoch_no <- proc_events(epoch_events = epoch_events,
                            event_type = data$events$event_type,
                            epoch_nos = data$events$epoch,
                            event_labels = data$events$event_label,
                            keep = keep)
  }
  if (is.numeric(epoch_no)) {
    if (keep == FALSE) {
      orig_epo_no <- unique(data$timings$epoch)
      epoch_no <- orig_epo_no[!orig_epo_no %in% epoch_no]
    }
    keep_rows <- data$timings$epoch %in% epoch_no
    data$signals <- data$signals[keep_rows, ]
    data$timings <- data$timings[keep_rows, ]
    data$events <- data$events[data$events$epoch %in% epoch_no, ]
  }
  if (df_out) {
    return(as.data.frame(data))
  }
  data
}

#' @describeIn select_epochs Selection of epochs from \code{eeg_ICA} objects.
#' @export

select_epochs.eeg_ICA <- function(data,
                                  epoch_events = NULL,
                                  epoch_no = NULL,
                                  keep = TRUE,
                                  df_out = FALSE,
                                  ...) {

  # First check if epoch_events has been passed; if it's numeric, select epochs
  # based on event_type. If it's a character vector, check if those labels exist
  # in the data.

  if (!is.null(epoch_events)) {
    epoch_no <- proc_events(epoch_events = epoch_events,
                            event_type = data$events$event_type,
                            epoch_nos = data$events$epoch,
                            event_labels = data$events$event_label,
                            keep = keep)
  }

  if (is.numeric(epoch_no)) {
    keep_rows <- data$timings$epoch %in% epoch_no
    if (keep == FALSE) {
      keep_rows <- !keep_rows
    }
    data$comp_activations <- data$comp_activations[keep_rows, ]
    data$timings <- data$timings[keep_rows, ]
    data$events <- data$events[data$events$epoch %in% epoch_no, ]
  }
  if (df_out) {
    return(as.data.frame(data))
  }
  data
}

#'@noRd
select_epochs.eeg_tfr <- function(data,
                                  epoch_events = NULL,
                                  epoch_no = NULL,
                                  keep = TRUE,
                                  df_out = FALSE,
                                  ...) {
  if ("epoch" %in% data$dimensions) {

  } else {
    stop("Data is averaged, so no epochs present.")
  }
  data
}

#' Select frequencies
#'
#' Select specific frequencies from \code{eeg_tfr} objects. Can be used to
#' selecting either single frequencies or anything within a range.
#'
#' @examples
#' demo_tfr <- compute_tfr(demo_epochs, foi = c(4, 30), n_freq = 10, n_cycles = 5)
#' select_freqs(demo_tfr, c(8, 12))
#' @param data An \code{eeg_tfr} object.
#' @param freq_range The range of frequencies to retain. Can be a scale or the
#'   lower and upper bounds. (e.g. c(5, 30))
#' @export
select_freqs <- function(data,
                         freq_range) {
  UseMethod("select_freqs", data)
}

#' @export
select_freqs.default <- function(data,
                                 freq_range) {

  warning(paste("select_freqs() does not handle objects of class", class(data),
                "and can currently only be used on eeg_tfr objects."))
}

#' @describeIn select_freqs Function for selecting specific frequencies from \code{eeg_tfr} objects.
#' @export
select_freqs.eeg_tfr <- function(data,
                                 freq_range) {

  freq_dim <- which(data$dimensions == "frequency")
  if (length(freq_range) == 2) {
    data_freqs <- as.numeric(dimnames(data$signals)[[freq_dim]])
    freqs <- data_freqs >= freq_range[[1]] &
      data_freqs <= freq_range[[2]]
    data$freq_info$freqs <- data$freq_info$freqs[freqs]
  } else if (length(freq_range) == 1) {

    closest_freq <- which.min(abs(data$freq_info$freqs - freq_range[1]))
    freqs <- closest_freq
    message(paste("Returning closest frequency, ", data$freq_info$freqs[freqs], "Hz"))
    data$freq_info$freqs <- data$freq_info$freqs[closest_freq]
  }

  data$signals <-
    abind::asub(data$signals,
                freqs,
                freq_dim,
                drop = FALSE)
  data
}

#' Internal function for processing epoch_events during selection
#'
#' Converts character strings into a vector of epoch numbers with matching labels.
#'
#' @keywords internal

proc_events <- function(epoch_events,
                        event_type,
                        epoch_nos,
                        event_labels,
                        keep
                        ) {

  if (is.numeric(epoch_events)) {
    keep_rows <- event_type %in% epoch_events
    if (!any(keep_rows)) {
      stop("Events not found.")
    }
    if (keep == FALSE) {
      keep_rows <- !keep_rows
    }
    epoch_no <- as.numeric(epoch_nos[keep_rows])
  } else if (is.character(epoch_events)) {
    check_ev <- label_check(epoch_events, event_labels)
    if (check_ev) {
      check_ev <- grepl(epoch_events, event_labels)
    } else {
      stop("Event label not found, check with list_events.")
    }
    keep_rows <- check_ev
    if (keep == FALSE) {
      keep_rows <- !keep_rows
    }
    epoch_no <- as.numeric(epoch_nos[keep_rows])
  }
}
