#' Select timerange
#'
#' Generic function for selecting specific time ranges from a given dataset.
#' Input can be a dataframe or an object of class \code{eeg_data}.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @param data Data from which to select
#' @param ... Further arguments passed to or from other methods.
#' @seealso \code{\link{select_times.default}},
#'   \code{\link{select_times.eeg_data}}, \code{\link{select_elecs}}
#'
#' @export

select_times <- function(data, ...) {
  UseMethod("select_times", data)
}

#' Select timerange
#'
#' Select timerange from a dataframe.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @param data An EEG dataset in a data frame. Must have a column named "time".
#' @param ... Arguments used with related methods
#' @param time_lim A character vector of two numbers indicating the time range to be selected e.g. c(min, max)
#' @importFrom dplyr filter
#' @return Data frame with only data from within the specified range.
#' @seealso \code{\link{select_times}}, \code{\link{select_times.eeg_data}}, \code{\link{select_elecs}}
#' @export

select_times.default <- function(data, time_lim = NULL, ...) {

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      warning("Must enter two timepoints when selecting a time range; using whole range.")
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

#' Select time range
#'
#' select a specifed time range from objects of class \code{eeg_data}.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @param data Must be an object of class \code{eeg_data}.
#' @param ... Arguments used with related methods
#' @param time_lim A character vector of two numbers indicating the time range
#'   to be selected e.g. c(min, max)
#' @param df_out Return a data frame rather than an \code{eeg_data} object.
#'   Defaults to FALSE (i.e. returns an \code{eeg_data} object)
#' @importFrom dplyr filter select
#' @export

select_times.eeg_data <- function(data, time_lim = NULL, df_out = FALSE, ...) {
  proc_data <- as.data.frame(data)

  if ("time" %in% colnames(proc_data)) {
    if (length(time_lim) == 1) {
      warning("Must enter two timepoints when selecting a time range; using whole range.")
    } else if (length(time_lim) == 2) {
      time_lim[1] <- proc_data$time[which.min(abs(proc_data$time - time_lim[1]))]
      time_lim[2] <- proc_data$time[which.min(abs(proc_data$time - time_lim[2]))]
      proc_data <- dplyr::filter(proc_data, time >= time_lim[1] & time <= time_lim[2])
    } else {
      warning("No time column found.")
    }
  }
  if (df_out) {
    return(proc_data)
  } else {
    data$signals <- dplyr::select(proc_data, -sample, -time)
    data$timings <- list(time = proc_data$time, sample = proc_data$sample)
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
#' @param electrode A character vector of electrode labels for selection or removal.
#' @param keep Defaults to TRUE. Set to false to *remove* the selected electrodes.
#'
#' @return Data frame with only data from the chosen electrodes
#'
#' @export
#'

select_elecs <- function(data, ...) {
  UseMethod("select_elecs", data)
}

#' Select electrodes from a given dataset.
#'
#' Checks for presence of electrode column, and if found, presence of selected electrodes.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @param data An EEG dataset.
#' @param ... Arguments used with related methods
#' @param electrode A character vector of electrode labels for selection or removal.
#' @param keep Defaults to TRUE. Set to false to *remove* the selected electrodes.
#'
#' @return Data frame with only data from the chosen electrodes
#'
#' @export
#'

select_elecs.default <- function(data,  electrode = NULL, keep = TRUE, ...) {

  if ("electrode" %in% colnames(data)) {
    if (all(electrode %in% data$electrode)) {
      if (keep) {
        data <- data[data$electrode %in% electrode, ]
      } else {
        data <- data[!data$electrode %in% electrode, ]
      }
    } else {
      warning(cat("Electrode(s) not found:", electrode[!electrode %in% data$electrode], ". Returning all data."))
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


#' Select electrodes from a given dataset.
#'
#' Checks for presence of electrode column, and if found, presence of selected electrodes.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @param data An \code{eeg_data} object.
#' @param ... Arguments used with related methods
#' @param electrode A character vector of electrode labels for selection or removal.
#' @param keep Defaults to TRUE. Set to false to *remove* the selected electrodes.
#' @param df_out Defaults to FALSE. Set to TRUE to return a dataframe rather than an \code{eeg_data} object.
#'
#' @return \code{eeg_data} object with selected electrodes removed/kept.
#'
#' @export
#'

select_elecs.eeg_data <-
  function(data,
           electrode,
           keep = TRUE,
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
