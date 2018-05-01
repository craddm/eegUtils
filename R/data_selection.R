#' Select timerange
#'
#' Generic function for selecting specific time ranges from a given dataset.
#' Input can be a dataframe, or an object of class \code{eeg_data} or
#' \code{eeg_epochs}. Note this finds the closest times to those specified, so
#' times out may not correspond exactly to requested times.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @param data Data from which to select
#' @param ... Further arguments passed to or from other methods.
#'
#' @export

select_times <- function(data, ...) {
  UseMethod("select_times", data)
}

#' @param time_lim A character vector of two numbers indicating the time range
#'   to be selected e.g. c(min, max)
#' @importFrom dplyr filter
#' @return Data frame with only data from within the specified range.
#' @export
#' @describeIn select_times Default select times function

select_times.default <- function(data, time_lim = NULL, ...) {

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      warning("Must enter two timepoints when selecting a time range;
              using whole range.")
    } else if (length(time_lim) == 2) {
      time_lim[1] <- data$time[which.min(abs(data$time - time_lim[1]))]
      time_lim[2] <- data$time[which.min(abs(data$time - time_lim[2]))]
      data <- dplyr::filter(data, time >= time_lim[1] & time <= time_lim[2])
    }
  } else {
    warning("No time column found.")
  }
  return(data)
}

#' @param df_out Returns a data frame rather than an object of the same type that was passed in
#' @importFrom dplyr filter select
#' @export
#'
#' @describeIn select_times Select times from an eeg_data object

select_times.eeg_data <- function(data, time_lim = NULL, df_out = FALSE, ...) {

  proc_data <- as.data.frame(data)
  proc_data <- select_times(proc_data, time_lim = time_lim)

  if (df_out) {
    return(proc_data)
  } else {
    data$events <- dplyr::filter(data$events, event_time >= time_lim[1],
                                 event_time <= time_lim[2])
    data$signals <- dplyr::select(proc_data, -sample, -time)
    data$timings <- tibble::tibble(time = proc_data$time, sample = proc_data$sample)
    if (!is.null(data$reference$ref_data)) {
      data$reference$ref_data <- data$reference$ref_data[data$timings$sample,]
    }
    return(data)
  }
}

#' @importFrom dplyr filter select
#' @export
#' @describeIn select_times Select times in \code{eeg_epoch} objects

select_times.eeg_epochs <- function(data, time_lim = NULL,
                                    df_out = FALSE, ...) {

  proc_data <- as.data.frame(data)
  proc_data <- select_times(proc_data, time_lim = time_lim)

  if (df_out) {
    return(proc_data)
  } else {
    data$events <- dplyr::filter(data$events, time >= time_lim[1],
                                 time <= time_lim[2])
    data$signals <- dplyr::select(proc_data, -sample, -time, -epoch)
    data$timings <- tibble::tibble(time = proc_data$time,
                         sample = proc_data$sample,
                         epoch = proc_data$epoch)

    if (!is.null(data$reference$ref_data)) {
      data$reference$ref_data <- data$reference$ref_data[data$timings$sample]
    }
    return(data)
  }
}

#' @export
#' @describeIn select_times Select times in \code{eeg_evoked} objects
select_times.eeg_evoked <- function(data, time_lim = NULL,
                                    df_out = FALSE, ...) {

  data$signals <- as.data.frame(data)
  data$signals <- select_times(data$signals, time_lim = time_lim)

  if (df_out) {
    return(data$signals)
  } else {
    return(data)
  }

}

#' Select electrodes from a given dataset.
#'
#' This is a generic function for selection of electrodes from an EEG dataset.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @param data An EEG dataset.
#' @param ... Arguments used with related methods
#'
#' @export
#'

select_elecs <- function(data, ...) {
  UseMethod("select_elecs", data)
}

#' @param electrode A character vector of electrode labels for selection or
#'   removal.
#' @param keep Defaults to TRUE. Set to false to *remove* the selected
#'   electrodes.
#'
#' @return Data frame with only data from the chosen electrodes
#'
#' @describeIn select_elecs Select electrodes from a generic data frame.
#' @export

select_elecs.default <- function(data,  electrode = NULL, keep = TRUE, ...) {

  if ("electrode" %in% colnames(data)) {
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
    if (all(electrode %in% colnames(data))) {
      if (keep) {
        data <- data[, colnames(data) %in% electrode, drop = FALSE]
      } else {
        data <- data[, !colnames(data) %in% electrode, drop = FALSE]
      }
    }
  }
  return(data)
}


#' @param df_out Defaults to FALSE. Set to TRUE to return a dataframe rather
#'   than an \code{eeg_data} object.
#'
#' @return \code{eeg_data} object with selected electrodes removed/kept.
#'
#' @export
#' @describeIn select_elecs Select electrodes from a \code{eeg_data} object.

select_elecs.eeg_data <- function(data, electrode, keep = TRUE,
                                  df_out = FALSE, ...) {

  if (all(electrode %in% colnames(data$signals))) {
    if (keep) {
      data$signals <- data$signals[colnames(data$signals) %in% electrode]
    } else {
      data$signals <- data$signals[!colnames(data$signals) %in% electrode]
    }
  } else {
    cat("Electrode(s) not found:",
        electrode[!electrode %in% colnames(data$signals)],
        ". Returning all data.")
    warning()
  }
  if (df_out) {
    return(as.data.frame(data))
  } else {
    return(data)
  }
}

#' Select epochs from eeg_data
#'
#' This is a generic function for selecting epochs from an epoched data set.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @param data \code{eeg_epochs} object from which to select epochs.
#' @param ... Parameters passed to specific methods
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
#'   be numeric or character vector. Will override any epoch_no input
#' @param epoch_no Select epochs by epoch number.
#' @param keep Defaults to TRUE, meaning select the specified epochs. Set to
#'   FALSE to remove specified epochs.
#' @param df_out Output a data.frame instead of an eeg_data object.
#' @describeIn select_epochs Selection of epochs from \code{eeg_epochs} objects.
#' @export

select_epochs.eeg_epochs <- function(data, epoch_events = NULL, epoch_no = NULL,
                                   keep = TRUE, df_out = FALSE, ...) {

  # First check if epoch_events has been passed; if it's numeric, select epochs
  # based on event_type. If it's a character vector, check if those labels exist
  # in the data.

  if (is.numeric(epoch_events)) {
    sel_rows <- data$events$event_type %in% epoch_events

    if (keep == FALSE) {
      sel_rows <- !sel_rows
    }

    epoch_no <- as.numeric(data$events$epoch[sel_rows])
  } else if (is.character(epoch_events)) {
    check_ev <- epoch_events %in% list_events(data)$event_label

    if (!all(check_ev)) {
      stop("Event label not found, check with list_events.")
    } else {
      epoch_events <- epoch_events[check_ev]
    }

    sel_rows <- data$events$event_label %in% epoch_events

    if (keep == FALSE) {
      sel_rows <- !sel_rows
    }

    epoch_no <- as.numeric(data$events$epoch[sel_rows])
  }

  if (is.numeric(epoch_no)) {
    sel_rows <- data$timings$epoch %in% epoch_no

    if (keep == FALSE) {
      sel_rows <- !sel_rows
    }

    data$signals <- data$signals[sel_rows, ]

    if (!is.null(data$reference)) {
      data$reference$ref_data <- data$reference$ref_data[sel_rows]
    }

    data$timings <- data$timings[sel_rows, ]
    data$events <- data$events[data$events$epoch %in% epoch_no, ]

  }

  if (df_out) {
    as.data.frame(data)
  } else {
    data
  }
}
