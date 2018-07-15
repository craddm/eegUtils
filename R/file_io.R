#' Function for reading raw data.
#'
#' Currently BDF/EDF and 32-bit .CNT files are supported. Filetype is determined
#' by the file extension.The \code{edfReader} package is used to load BDF/EDF files.
#' The function creates an \code{eeg_data} structure for subsequent use.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param file_name File to import. Should include file extension.
#' @param file_path Path to file name, if not included in filename.
#' @param chan_nos Channels to import. All channels are included by default.
#' @import edfReader
#' @import tools
#' @importFrom purrr map_df is_empty
#' @importFrom tibble tibble as_tibble
#' @export

import_raw <- function(file_name, file_path = NULL, chan_nos = NULL) {
  file_type <- tools::file_ext(file_name)

  if (file_type == "bdf" | file_type == "edf") {
    data <- edfReader::readEdfSignals(edfReader::readEdfHeader(file_name))

    #check for an annotations channel
    anno_chan <- which(vapply(data, function(x) isTRUE(x$isAnnotation),
                 FUN.VALUE = logical(1)))

    #remove annotations if present - could put in separate list...
    if (length(anno_chan) > 0) {
      data <- data[-anno_chan]
    }

    sigs <- purrr::map_df(data, "signal")
    srate <- data[[1]]$sRate

    if ("Status" %in% names(sigs)) {
      events <- sigs$Status %% (256)
    } else {
      events <- integer(0)
    }

    if (purrr::is_empty(events)) {
      warning("No status channel. Retaining all channels.")
      chan_nos <- 1:ncol(sigs)
    } else {
      #all chans except Status
      chan_nos <- (1:ncol(sigs))[-which(names(sigs) == "Status")]
    }

    timings <- tibble::tibble(sample = 1:nrow(sigs))
    timings$time <- (timings$sample - 1) / srate

    if (is.null(chan_nos)) {
      chan_nos <- 1:(ncol(sigs) - 1)
    }

    sigs <- tibble::as_tibble(sigs[, chan_nos])
    events_diff <- diff(events)
    event_table <- tibble::tibble(event_onset = which(events_diff > 0) + 1,
                                  event_time = which(events_diff > 0) / srate,
                              event_type = events[which(events_diff > 0) + 1])
    data <- eeg_data(data = sigs, srate = srate,
                     events = event_table, timings = timings,
                     continuous = TRUE)
  } else if (file_type == "cnt") {
    data <- import_cnt(file_name)
    sigs <- tibble::tibble(t(data$chan_data))
    names(sigs) <- data$chan_info$electrode
    srate <- data$head_info$samp_rate
    timings <- tibble::tibble(sample = 1:dim(sigs)[[1]])
    timings$time <- (timings$sample - 1) / srate
    event_table <- tibble::tibble(event_onset = data$event_list$offset + 1,
                                  event_time = (data$event_list$offset + 1) / srate,
                                  event_type = data$event_list$event_type)
    data <- eeg_data(data = sigs, srate = srate,
                     chan_info = data$chan_info,
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
#' Beta version of function to import Neuroscan .CNT files. Only intended for
#' import of 32-bit files.
#'
#' @param file_name Name of .CNT file to be loaded.
#' @importFrom tibble tibble
#' @noRd

import_cnt <- function(file_name) {

  cnt_file <- file(file_name, "rb")

  # Read in meta-info - number of channels, location of event table, sampling
  # rate...

  pos <- seek(cnt_file, 12)
  next_file <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
  pos <- seek(cnt_file, 353)
  n_events <- readBin(cnt_file, integer(), n = 1, endian = "little")
  pos <- seek(cnt_file, 370)
  n_channels <- readBin(cnt_file, integer(), n = 1, size = 2,
                        signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 376)
  samp_rate <-  readBin(cnt_file, integer(), n = 1, size = 2,
                        signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 864)
  n_samples <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
  pos <- seek(cnt_file, 886)
  event_table_pos <- readBin(cnt_file, integer(), size = 4,
                             n = 1, endian = "little") # event table
  pos <- seek(cnt_file, 900)

  data_info <- data.frame(n_events,
                          n_channels,
                          samp_rate,
                          event_table_pos)

  chan_df <- data.frame(electrode = character(n_channels),
                        chan_no = numeric(n_channels),
                        x = numeric(n_channels),
                        y = numeric(n_channels),
                        baseline = numeric(n_channels),
                        sens = numeric(n_channels),
                        cal = numeric(n_channels),
                        stringsAsFactors = FALSE)

  # Read channel names and locations

  for (i in 1:n_channels) {
    chan_start <- seek(cnt_file)
    chan_df$electrode[i] <- readBin(cnt_file, character(),
                                    n = 1, endian = "little")
    chan_df$chan_no[i] <- i
    pos <- seek(cnt_file, chan_start + 19)
    chan_df$x[i] <- readBin(cnt_file, double(), size = 4,
                            n = 1, endian = "little") # x coord
    chan_df$y[i] <- readBin(cnt_file, double(), size = 4,
                            n = 1, endian = "little") # y coord

    pos <- seek(cnt_file, chan_start + 47)
    chan_df$baseline[i] <- readBin(cnt_file, integer(), size = 1,
                                   n = 1, endian = "little")
    pos <- seek(cnt_file, chan_start + 59)
    chan_df$sens[i] <- readBin(cnt_file, double(), size = 4, n = 1,
                               endian = "little")
    pos <- seek(cnt_file, chan_start + 71)
    chan_df$cal[i] <- readBin(cnt_file, double(), size = 4, n = 1,
                              endian = "little")
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

  ev_list <- data.frame(event_type = integer(n_events),
                        keyboard = character(n_events),
                        keypad_accept = integer(n_events),
                        accept_evl = integer(n_events),
                        offset = integer(n_events),
                        type = integer(n_events),
                        code = integer(n_events),
                        latency = numeric(n_events),
                        epochevent = integer(n_events),
                        accept = integer(n_events),
                        accuracy = integer(n_events),
                        stringsAsFactors = FALSE
                        )

  for (i in 1:n_events) {
    ev_list$event_type[i] <- readBin(cnt_file, integer(), size = 2,
                                     n = 1, endian = "little")
    ev_list$keyboard[i] <- readBin(cnt_file, integer(), size = 1,
                                   n = 1, endian = "little")
    temp <- readBin(cnt_file, integer(), size = 1, n = 1, signed = FALSE,
                    endian = "little")
    ev_list$keypad_accept[i] <- bitwAnd(15, temp)
    ev_list$accept_evl[i] <- bitwShiftR(temp, 4)
    ev_list$offset[i] <- readBin(cnt_file, integer(), size = 4,
                                 n = 1, endian = "little")
    ev_list$type[i] <- readBin(cnt_file, integer(), size = 2,
                               n = 1, endian = "little")
    ev_list$code[i] <- readBin(cnt_file, integer(), size = 2,
                               n = 1, endian = "little")
    ev_list$latency[i] <- readBin(cnt_file, double(), size = 4,
                                  n = 1, endian = "little")
    ev_list$epochevent[i] <- readBin(cnt_file, integer(), size = 1,
                                     n = 1, endian = "little")
    ev_list$accept[i] <- readBin(cnt_file, integer(), size = 1,
                                 n = 1, endian = "little")
    ev_list$accuracy[i] <- readBin(cnt_file, integer(), size = 1,
                                   n = 1, endian = "little")
  }

  ev_list$offset <- (ev_list$offset - beg_data) / (4 * n_channels) + 1

  close(cnt_file)
  chan_df <- chan_df[, names(chan_df) %in% c("electrode", "chan_no")]
  out <- list(chan_info = chan_df, head_info = data_info,
              chan_data = chan_data, event_list = ev_list)
}

#' Load EEGLAB .set files
#'
#' DEPRECATED - use \code{import_set}
#'
#' EEGLAB .set files are standard Matlab .mat files, but EEGLAB can be set to
#' export either v6.5 or v7.3 format files. Only v6.5 files can be read with
#' this function. v7.3 files (which use HDF5 format) are not currently
#' supported, as they cannot be fully read with existing tools.
#'
#'
#' @param file_name Filename (and path if not in present working directory)
#' @param df_out Defaults to FALSE - outputs an object of class eeg_data. Set to
#'   TRUE for a normal data frame.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @importFrom R.matlab readMat
#' @importFrom dplyr group_by mutate rename
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr is_empty
#' @export

load_set <- function(file_name, df_out = FALSE) {
  .Deprecated("import_set",
              msg = c("load_set is deprecated and has been renamed to import_set. ",
              "It will be removed in a forthcoming release."))
  file_out <- import_set(file_name, df_out = FALSE)
  file_out
}

#' Load EEGLAB .set files
#'
#' EEGLAB .set files are standard Matlab .mat files, but EEGLAB can be set to
#' export either v6.5 or v7.3 format files. Only v6.5 files can be read with
#' this function. v7.3 files (which use HDF5 format) are not currently
#' supported, as they cannot be fully read with existing tools.
#'
#' @param file_name Filename (and path if not in present working directory)
#' @param df_out Defaults to FALSE - outputs an object of class eeg_data. Set to
#'   TRUE for a normal data frame.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @importFrom R.matlab readMat
#' @importFrom dplyr group_by mutate rename
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr is_empty
#' @export

import_set <- function(file_name, df_out = FALSE) {

  temp_dat <- R.matlab::readMat(file_name)
  var_names <- dimnames(temp_dat$EEG)[[1]]

  n_chans <- temp_dat$EEG[[which(var_names == "nbchan")]]
  n_trials <- temp_dat$EEG[[which(var_names == "trials")]]
  times <- temp_dat$EEG[[which(var_names == "times")]]

  chan_info <- temp_dat$EEG[[which(var_names == "chanlocs")]]
  col_names <- dimnames(chan_info)[1]
  size_chans <- dim(chan_info)
  chan_info <- lapply(chan_info,
                      function(x) ifelse(purrr::is_empty(x), NA, x))
  dim(chan_info) <- size_chans
  dimnames(chan_info) <- col_names
  chan_info <- parse_chaninfo(chan_info)

  # check if the data is stored in the set or in a separate .fdt
  if (is.character(temp_dat$EEG[[which(var_names == "data")]])) {
    message("loading from .fdt")
    fdt_file <- paste0(tools::file_path_sans_ext(file_name), ".fdt")
    fdt_file <- file(fdt_file, "rb")

    # read in data from .fdt
    # do this in chunks to avoid memory errors for large files...?
    signals <- readBin(fdt_file,
                       "double",
                       n = n_chans * n_trials * length(times),
                       size = 4,
                       endian = "little")
    close(fdt_file)

    dim(signals) <- c(n_chans, length(times) * max(n_trials, 1))
    times <- rep(times, max(n_trials, 1))

    if (n_trials == 1) {
      continuous <- TRUE
    } else {
      continuous <- FALSE
    }

  } else {

    # if the data is in the .set file, load it here instead of above
    signals <- temp_dat$EEG[[which(dimnames(temp_dat$EEG)[[1]] == "data")]]
    dim_signals <- dim(signals)

    if (length(dim_signals) == 3) {
      dim(signals) <- c(dim_signals[1], dim_signals[2] * dim_signals[3])
      times <- rep(times, n_trials)
      continuous <- FALSE
    } else {
      continuous <- TRUE
    }
  }

  signals <- data.frame(cbind(t(signals), times))
  srate <- c(temp_dat$EEG[[which(var_names == "srate")]])
  names(signals) <- c(unique(chan_info$electrode), "time")
  signals <- dplyr::group_by(signals, time)
  signals <- dplyr::mutate(signals, epoch = 1:n())
  signals <- dplyr::ungroup(signals)

  event_info <- temp_dat$EEG[[which(var_names == "event")]]

  event_table <- tibble::as_tibble(t(matrix(as.integer(event_info),
                               nrow = dim(event_info)[1],
                               ncol = dim(event_info)[3])))

  names(event_table) <- unlist(dimnames(event_info)[1])
  event_table$event_time <- (event_table$latency - 1) / srate
  event_table <- event_table[c("latency", "event_time", "type", "epoch")]
  event_table <- dplyr::rename(event_table, event_type = "type",
                               event_onset = "latency")
  event_table$time <- NA

  if (df_out) {
    return(signals)
  } else {
    signals$time <- signals$time / 1000
    # convert to seconds - eeglab uses milliseconds
    timings <- tibble::tibble(time = signals$time,
                              epoch = signals$epoch,
                              sample = 1:length(signals$time))
    event_table$time <- timings[which(timings$sample %in% event_table$event_onset, arr.ind = TRUE), ]$time
    out_data <- eeg_data(signals[, 1:n_chans],
                         srate = srate,
                         timings = timings,
                         continuous = continuous,
                         chan_info = chan_info,
                         events = event_table)
    if (!continuous) {
      class(out_data) <- c("eeg_epochs", "eeg_data")
    }
    out_data
  }
}


#' Parse channel info from an EEGLAB set file
#'
#' Internal function to convert EEGLAB chan_info to eegUtils style
#'
#' @param chan_info Channel info list from an EEGLAB set file
#' @noRd
parse_chaninfo <- function(chan_info) {
  chan_info <- tibble::as_tibble(t(as.data.frame(chan_info)))
  chan_info <- tibble::as_tibble(data.frame(lapply(chan_info,
                                                   unlist),
                                            stringsAsFactors = FALSE))
  names(chan_info) <- c("electrode",
                        "angle",
                        "radius",
                        "cart_x",
                        "cart_y",
                        "cart_z",
                        "sph_theta",
                        "sph_phi",
                        "sph_radius",
                        "type",
                        "ref",
                        "urchan")
  #pol_coords <- cart_to_pol(chan_info$cart_x, chan_info$cart_y)
  #chan_info <- tibble::as_tibble(c(chan_info, pol_coords))
  chan_info <- chan_info[c("electrode",
                           "cart_x",
                           "cart_y",
                           "cart_z",
                           "sph_radius",
                           "sph_phi",
                           "sph_theta",
                           "angle",
                           "radius")]
  # EEGLAB co-ordinates are rotated 90 degrees compared to our coordinates,
  # and left-right flipped
  chan_info <- rotate_sph(chan_info, -90)
  chan_info <- flip_x(chan_info)
  chan_info
  }
