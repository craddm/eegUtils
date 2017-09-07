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
      data <- dplyr::filter(data, time >= time_lim[1] & time <= time_lim[2])
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
