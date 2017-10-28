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


#' Import Neuroscan .CNT file
#'
#'
#' @param file_name Name of .CNT file to be loaded.

import_cnt <- function(file_name) {
  cnt_file <- file(file_name, "rb")
  pos <- seek(cnt_file, 353)
  n_events <- readBin(cnt_file, integer(), n = 1, endian = "little")
  pos <- seek(cnt_file, 370)
  n_channels <- readBin(cnt_file, integer(), n = 1, size = 2, signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 376)
  samp_rate <-  readBin(cnt_file, integer(), n = 1, size = 2, signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 886)
  event_table_pos <- readBin(conn, integer(), size= 4,  n = 1, endian = "little") # event table
  pos <- seek(cnt_file, 900)
  chan_details <- vector("list", n_channels)
  chan_df <- tibble::tibble(chan_name = character(n_channels),
                            chan_no = numeric(n_channels),
                            x = numeric(n_channels),
                            y = numeric(n_channels)
                            )

  for (i in 1:n_channels) {
    chan_start <- seek(cnt_file)
    chan_name <- readBin(cnt_file, character(), n = 1, endian = "little")
    pos <- seek(cnt_file, chan_start + 19)
    xcoord <- readBin(conn, double(), size = 4,  n = 1, endian = "little") # x coord
    ycoord <- readBin(conn, double(), size = 4,  n = 1, endian = "little") # y coord

    chan_df$chan_name[i] <- chan_name
    chan_df$chan_no[i] <- i
    chan_df$x[i] <- xcoord
    chan_df$y[i] <- ycoord
    pos <- seek(cnt_file, (900 + i * 75))
  }
  beg_data <- seek(cnt_file) # beginning of actual data
  n_samples <- event_table_pos - (900 + 75 * n_channels) / (2 * n_channels)

  close(cnt_file)

  chan_df
}
