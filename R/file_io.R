#' Function for reading raw data.
#'
#' Currently only BDF or EDF files are supported. The \code{edfReader} package reads the data in, and then the function creates an eeg_data structure for subsequent use.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @param file_name File to import.
#' @param file_path Path to file name, if not included in filename.
#' @param chan_nos Channels to import. All channels are included by default.
#' @import edfReader
#' @import tools
#' @importFrom purrr map_df
#' @importFrom tibble tibble as_tibble
#' @export

import_raw <- function(file_name, file_path = NULL, chan_nos = NULL) {
  file_type <- tools::file_ext(file_name)
  if (file_type == "bdf" | file_type == "edf") {
    data <- edfReader::readEdfSignals(edfReader::readEdfHeader(file_name))
    sigs <- purrr::map_df(data, "signal")
    srate <- data[[1]]$sRate
    events <- sigs$Status %% (256)
    timings <- tibble::tibble(sample = 1:dim(sigs)[[1]])
    timings$time <- (timings$sample - 1) / srate
    if (is.null(chan_nos)) {
      chan_nos <- 1:(dim(sigs)[[2]] - 1)
    }
    sigs <- tibble::as_tibble(sigs[,chan_nos])
    events_diff <- diff(events)
    event_table <- tibble::tibble(event_onset = which(events_diff > 0) + 1,
                              event_type = events[which(events_diff > 0) + 1])
    data <- eeg_data(data = sigs, srate = srate,
                     events = event_table, timings = timings,
                     continuous = TRUE)
  } else {
    warning("Unsupported filetype")
    return()
  }
  return(data)
}
