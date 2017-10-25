#' Get standard electrode locations
#'
#' Joins standard electrode locations to EEG data from eegUtils internal data.
#'
#' @param data An EEG dataset.
#' @param electrode The column name containing electrode names in data.
#'   (Defaults to "electrode").
#' @param drop Should electrodes in \code{data} for which default locations are
#'   not available be dropped? (Defaults to FALSE).
#' @param plot Plot obtained electrode locations.
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble is.tibble
#' @return A tibble (or data.frame), or ggplot2 object if \code{plot = TRUE}.
#' @export

electrode_locations <- function(data,
                                electrode = "electrode",
                                drop = FALSE,
                                plot = FALSE,
                                montage = NULL) {

  if (!is.null(montage)) {
    if (montage == "biosemi64alpha") {
      electrodeLocs[1:64, electrode] <- c(paste0("A", 1:32), paste0("B", 1:32))
    }
  }

  data[, electrode] <- toupper(data[[electrode]])
  electrodeLocs[, electrode] <- toupper(electrodeLocs[[electrode]])

  if (tibble::is.tibble(data)) {
    elecs <- dplyr::pull(unique(data[,electrode])) %in% dplyr::pull(electrodeLocs[,electrode])
    if (!all(elecs)) {
      message("Electrodes not found: ", paste(dplyr::pull(unique(data[,electrode]))[!elecs], sep = ","))
      } else if (!any(elecs)) {
        stop("No matching electrodes found.")
      }
  } else {
    elecs <- unique(data[,electrode]) %in% dplyr::pull(electrodeLocs[,electrode])
    if (!all(elecs)) {
      message("Electrodes not found: ", paste(unique(data[,electrode])[!elecs], sep = ","))
    } else if (!any(elecs)) {
      stop("No matching electrodes found.")
    }

  }

  if (drop) {
    data <- dplyr::inner_join(data, electrodeLocs, by = electrode)
  } else {
    data <- dplyr::left_join(data, electrodeLocs, by = electrode)
  }

  if (plot) {
    plotdata <- dplyr::distinct(data, x, y, electrode)
    p <- ggplot2::ggplot(plotdata, aes(x, y)) +
      geom_label(aes(label = electrode))
    return(p)
  } else {
    return(data)
  }
}

#' Function to create an S3 object of class "eeg_data".
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Raw data - signals from electrodes/channels.
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param continuous Whether the data is continuous or epoched.
#' @param reference Reference channel information, including names of reference channels, excluded channels etc.
#' @param events Event table
#' @export

eeg_data <- function(data,
                     srate,
                     events = NULL,
                     chan_info = NULL,
                     timings = NULL,
                     continuous = NULL,
                     reference = NULL) {
  if (srate < 1) {
    stop("Sampling rate must be above 0")
  }
  value <- list(signals = data,
                srate = srate,
                events = events,
                chan_info = chan_info,
                timings = timings,
                continuous = continuous,
                reference = reference
                )
  class(value) <- "eeg_data"
  value
}

#' Check if object is of class "eeg_data".
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object to check.
#'

is.eeg_data <- function(x) inherits(x, "eeg_data")

#' Convert eeg_data to data.frame
#'
#' Convert an object of class \code{eeg_data} into a standard data.frame / tibble
#'
#' @param data Object of class \code{eeg_data}
#' @param ... arguments for other as.data.frame commands
#' @param long Convert to long format. Defaults to FALSE
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_data <- function (data, ..., long = FALSE) {
  df <- data.frame(data$signals, data$timings)
  if (long) {
    if (data$continuous) {
      df <- tidyr::gather(df, electrode, amplitude, -time, -sample)
    } else {
      df <- tidyr::gather(df, electrode, amplitude, -time, -sample, -epoch)
    }
  }
  return(df)
}

#' Switch from wide to long format.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to convert
#' @importFrom tidyr gather

switch_format <- function(data) {
  data <- tidyr::gather(data, electrode, amplitude, -time)
  if (is.eeg_data(data)) {
  }
}
