#' Get standard electrode locations
#'
#' Joins standard electrode locations to EEG data from eegUtils internal data.
#'
#' @param data An EEG dataset.
#' @param electrode The column name containing electrode names in data.
#' (Defaults to "electrode").
#' @param drop Should electrodes in \code{data} for which default locations
#' are not available be dropped? (Defaults to FALSE).
#' @param plot Plot obtained electrode locations.
#'
#' @import dplyr
#' @return A tibble (or data.frame), or ggplot2 object if \code{plot = TRUE}.
#' @export

electrode_locations <- function(data,
                                electrode = "electrode",
                                drop = FALSE,
                                plot = FALSE) {
  data[, electrode] <- toupper(data[[electrode]])
  electrodeLocs[, electrode] <- toupper(electrodeLocs[[electrode]])

  if (drop) {
    data <- inner_join(data, electrodeLocs, by = electrode)
  } else {
    data <- left_join(data, electrodeLocs, by = electrode)
  }

  if (plot) {
    plotdata <- distinct(data, x, y, electrode)
    p <- ggplot(plotdata, aes(x, y)) +
      geom_label(aes(label = electrode))
    return(p)
  } else {
    return(data)
  }

}

#' Select timepoints from a given dataset
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#'
#' @param data An EEG dataset.
#' @param time_lim A character vector of two numbers indicating the time range to be selected e.g. c(min, max)
#' @import dplyr
#' @return Data frame with only data from within the specified range.
#'
#' @export

select_times <- function(data, time_lim) {

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      warning("Must enter two timepoints when selecting a time range; using whole range.")
    } else if (length(time_lim) == 2) {
      time_lim[1] <- data$time[which.min(abs(data$time - time_lim[1]))]
      time_lim[2] <- data$time[which.min(abs(data$time - time_lim[2]))]
      data <- filter(data, time >= time_lim[1] & time <= time_lim[2])
    } else {
      warning("No time column found.")
    }
  }
  return(data)
}


#' Select electrodes from a given dataset.
#'
#' Checks for presence of electrode column, and if found, presence of selected electrodes.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#'
#' @param data An EEG dataset.
#' @param electrode A character vector of electrode labels for selection or removal.
#' @param keep Defaults to TRUE. Set to false to *remove* the selected electrodes.
#'
#' @return Data frame with only data from the chosen electrodes
#'
#' @export
#'

select_elecs <- function(data, electrode, keep = TRUE) {
  if ("electrode" %in% colnames(data)) {
    if (all(electrode %in% data$electrode)) {
      if (keep) {
        data <- data[data$electrode %in% electrode, ]
      } else {
        data <- data[!data$electrode %in% electrode, ]
      }
      } else {
        cat("Electrode(s) not found:", electrode[!electrode %in% data$electrode], ". Returning all data.")
        warning()
      }
    } else {
    warning("No electrode column found.")
    }
  return(data)
}

#' Function to create an S3 object of class "eeg_data".
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data Raw data.
#' @param srate Sampling rate.
#' @param chan_labels
#'

eeg_data <- function(data, srate, events = NULL, chan_labels = NULL) {
  if (srate < 1) {
    stop("Sampling rate must be above 0")
  }
  value <- list(signals = data, srate = srate, events = events)
  class(value) <- "eeg_data"
  value
}

#' Check if object is of class "eeg_data".
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param x Object to check.
#'

is.eeg_data <- function(x) inherits(x, "eeg_data")


#' Switch from wide to long format.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data Data to convert
#' @import tidyr

switch_format <- function(x) {
  x <- gather(x, chan_label, amplitude, -time)
  if (is.eeg_data(x)) {
  }
}


