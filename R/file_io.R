#' Function for reading raw data.
#'
#' Currently BDF/EDF and 32-bit .CNT files are supported. Filetype is determined
#' by the file extension.The \code{edfReader} package is used to load BDF files.
#' The function creates an eeg_data structure for subsequent use.
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
  } else if (file_type == "cnt") {
    data <- import_cnt(file_name)
    sigs <- tibble::as_tibble(t(data$chan_data))
    names(sigs) <- data$chan_info$chan_name
    srate <- data$head_info$samp_rate
    timings <- tibble::tibble(sample = 1:dim(sigs)[[1]])
    timings$time <- (timings$sample - 1) / srate
    event_table <- tibble::tibble(event_onset = data$event_list$offset,
                                  event_type = data$event_list$event_type)
    data <- eeg_data(data = sigs, srate = srate, chan_info = data$chan_info[1:4],
                     events = event_table, timings = timings,
                     continuous = TRUE)
    } else{
    warning("Unsupported filetype")
    return()
  }
  return(data)
}


#' Import Neuroscan .CNT file
#'
#' Beta version of function to import Neuroscan .CNT files. Only intended for import of 32-bit files.
#'
#' @param file_name Name of .CNT file to be loaded.
#' @importFrom tibble tibble

import_cnt <- function(file_name) {

  cnt_file <- file(file_name, "rb")

  # Read in meta-info - number of channels, location of event table, sampling rate...

  pos <- seek(cnt_file, 12)
  next_file <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
  pos <- seek(cnt_file, 353)
  n_events <- readBin(cnt_file, integer(), n = 1, endian = "little")
  pos <- seek(cnt_file, 370)
  n_channels <- readBin(cnt_file, integer(), n = 1, size = 2, signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 376)
  samp_rate <-  readBin(cnt_file, integer(), n = 1, size = 2, signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 864)
  n_samples <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
  pos <- seek(cnt_file, 886)
  event_table_pos <- readBin(cnt_file, integer(), size= 4,  n = 1, endian = "little") # event table
  pos <- seek(cnt_file, 900)

  data_info <- tibble::tibble(n_events,
                              n_channels,
                              samp_rate,
                              event_table_pos)

  chan_df <- tibble::tibble(chan_name = character(n_channels),
                            chan_no = numeric(n_channels),
                            x = numeric(n_channels),
                            y = numeric(n_channels)
  )

  # Read channel names and locations

  for (i in 1:n_channels) {
    chan_start <- seek(cnt_file)
    chan_df$chan_name[i] <- readBin(cnt_file, character(), n = 1, endian = "little")
    chan_df$chan_no[i] <- i
    pos <- seek(cnt_file, chan_start + 19)
    chan_df$x[i] <- readBin(cnt_file, double(), size = 4,  n = 1, endian = "little") # x coord
    chan_df$y[i] <- readBin(cnt_file, double(), size = 4,  n = 1, endian = "little") # y coord

    pos <- seek(cnt_file, chan_start + 47)
    chan_df$baseline[i] <- readBin(cnt_file, integer(), size = 1, n = 1, endian = "little")
    pos <- seek(cnt_file, chan_start + 59)
    chan_df$sens[i] <- readBin(cnt_file, double(), size = 4, n = 1, endian = "little")
    pos <- seek(cnt_file, chan_start + 71)
    chan_df$cal[i] <- readBin(cnt_file, double(), size = 4, n = 1, endian = "little")
    pos <- seek(cnt_file, (900 + i * 75))
  }

  beg_data <- seek(cnt_file) # beginning of actual data
  real_n_samples <- event_table_pos - (900 + 75 * n_channels) / (2 * n_channels)
  frames <- floor((event_table_pos - beg_data) / n_channels / 4)

  chan_data <- matrix(readBin(cnt_file,
                              integer(),
                              size = 4,
                              n = n_channels * frames,
                              endian = "little"),
                      nrow = n_channels, ncol = frames)

  # rescale chan_data to microvolts
  mf <- chan_df$sens * (chan_df$cal / 204.8)
  chan_data <- (chan_data - chan_df$baseline) * mf

  # Read event table

  pos <- seek(cnt_file, event_table_pos)
  teeg <- readBin(cnt_file, integer(), size = 1, n = 1, endian = "little")
  tsize <- readBin(cnt_file, integer(), n = 1, endian = "little")
  toffset <- readBin(cnt_file, integer(), n = 1, endian = "little")
  ev_table_start <- seek(cnt_file)

  ev_list <- tibble::tibble(event_type = integer(n_events),
                            keyboard = character(n_events),
                            keypad_accept = integer(n_events),
                            accept_evl = integer(n_events),
                            offset = integer(n_events),
                            type = integer(n_events),
                            code = integer(n_events),
                            latency = numeric(n_events),
                            epochevent = integer(n_events),
                            accept = integer(n_events),
                            accuracy = integer(n_events)
                            )

  for (i in 1:n_events) {
    ev_list$event_type[i] <- readBin(cnt_file, integer(), size = 2, n = 1, endian = "little")
    ev_list$keyboard[i] <- readBin(cnt_file, integer(), size = 1, n = 1, endian = "little")
    temp <- readBin(cnt_file, integer(), size = 1, n = 1, signed = FALSE, endian = "little")
    ev_list$keypad_accept[i] <- bitwAnd(15, temp)
    ev_list$accept_evl[i] <- bitwShiftR(temp, 4)
    ev_list$offset[i] <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
    ev_list$type[i] <- readBin(cnt_file, integer(), size = 2, n = 1, endian = "little")
    ev_list$code[i] <- readBin(cnt_file, integer(), size = 2, n = 1, endian = "little")
    ev_list$latency[i] <- readBin(cnt_file, double(), size = 4, n = 1, endian = "little")
    ev_list$epochevent[i] <- readBin(cnt_file, integer(), size = 1,  n = 1, endian = "little")
    ev_list$accept[i] <- readBin(cnt_file, integer(), size = 1, n = 1, endian = "little")
    ev_list$accuracy[i] <- readBin(cnt_file, integer(),size = 1,  n = 1, endian = "little")
  }

  ev_list$offset <- (ev_list$offset - beg_data) / (4 * n_channels)

  close(cnt_file)

  #list(chan_df, data_info, chan_data)
  out <- list(chan_info = chan_df, head_info = data_info, chan_data = chan_data, event_list = ev_list)
}


#' Load EEGLAB .set files
#'
#' EEGLAB .set files are standard Matlab .mat files, but EEGLAB can be set to
#' export either v6.5 or v7.3 format files. Only v6.5 files can be read with
#' this function. v7.3 files (which use HDF5 format) are not currently
#' supported, as they cannot be read well with existing tools.
#'
#' @param file_name Filename (and path if not in present working directory)
#' @param df_out Defaults to FALSE - outputs an object of class eeg_data. Set to TRUE for a normal data frame.
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @import R.matlab
#' @importFrom dplyr group_by mutate


load_set <- function (file_name, df_out = FALSE) {
  file_type <- tools::file_ext(file_name)
  temp_dat <- R.matlab::readMat(file_name)
  n_chans <- temp_dat$EEG[[9]]
  n_trials <- temp_dat$EEG[[10]]
  times <- temp_dat$EEG[[15]]
  chan_info <- as.data.frame(temp_dat$EEG[[22]])
  if (is.character(temp_dat$EEG[[16]])) {
    message("loading from .fdt")
    fdt_file <- paste0(tools::file_path_sans_ext(file_name), '.fdt')
    fdt_file <- file(fdt_file, "rb")
    signals <-
      readBin(
        fdt_file,
        double(),
        n = n_chans * n_trials * length(times),
        size = 4,
        endian = "little"
      )
    close(fdt_file)
    dim(signals) <- c(n_chans, length(times) * max(n_trials, 1))
    times <- rep(times, max(n_trials, 1))
    if (n_trials == 1) {
      continuous <- TRUE
    } else {
      continuous <- FALSE
    }
  } else {
    signals <- temp_dat$EEG[[16]]
    dim_signals <- dim(signals)
    if (length(dim_signals) == 3) {
      dim(signals) <- c(dim_signals[1], dim_signals[2] * dim_signals[3])
      times <- rep(times, n_trials)
      continuous <- FALSE
    } else {
      continuous <- TRUE
    }
  }

  final_dat <- data.frame(cbind(t(signals), times))
  srate <- temp_dat$EEG[[12]]
  names(final_dat) <- c(unlist(chan_info[1,]), "time")
  final_dat <- dplyr::group_by(final_dat, time)
  final_dat <- dplyr::mutate(final_dat, epoch = 1:n())
  if (df_out) {
    return(final_dat)
  } else {
    final_dat$time <- final_dat$time / 1000 # convert to seconds
    timings <- data.frame(time = final_dat$time, epoch = final_dat$epoch)
    eeg_data(final_dat[, 1:n_chans], srate = srate, timings = timings, continuous = continuous)
  }
}
