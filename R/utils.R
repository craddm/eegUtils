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
