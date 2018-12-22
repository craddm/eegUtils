#' Function for reading raw data.
#'
#' Currently BDF/EDF, 32-bit .CNT, and Brain Vision Analyzer files are
#' supported. Filetype is determined by the file extension.The \code{edfReader}
#' package is used to load BDF/EDF files. The function creates an
#' \code{eeg_data} structure for subsequent use.
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
    anno_chan <- which(vapply(data,
                              function(x) isTRUE(x$isAnnotation),
                 FUN.VALUE = logical(1)))

    #remove annotations if present - could put in separate list...
    if (length(anno_chan) > 0) {
      data <- data[-anno_chan]
      message("Annotations are currently discarded. File an issue if you'd like
 this to change.")
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

    epochs <- data.frame(epoch = 1,
                         recording = basename(tools::file_path_sans_ext(file_name)))
    data <- eeg_data(data = sigs,
                     srate = srate,
                     events = event_table,
                     timings = timings,
                     continuous = TRUE,
                     epochs = epochs)
    data
  } else if (file_type == "cnt") {
    data <- import_cnt(file_name)
    sigs <- tibble::tibble(t(data$chan_data))
    names(sigs) <- data$chan_info$electrode
    srate <- data$head_info$samp_rate
    timings <- tibble::tibble(sample = 1:dim(sigs)[[1]])
    timings$time <- (timings$sample - 1) / srate
    event_table <-
      tibble::tibble(event_onset = data$event_list$offset + 1,
                     event_time = (data$event_list$offset + 1) / srate,
                     event_type = data$event_list$event_type)
    epochs <- data.frame(epoch = 1,
                         recording = basename(tools::file_path_sans_ext(file_name)))
    data <- eeg_data(data = sigs,
                     srate = srate,
                     chan_info = data$chan_info,
                     events = event_table,
                     timings = timings,
                     continuous = TRUE,
                     epochs = epochs)
    } else if (file_type == "vhdr") {
      data <- import_vhdr(file_name)
    } else {
      stop("Unsupported filetype")
    }
  data
}


#' Import Neuroscan .CNT file
#'
#' Beta version of function to import Neuroscan .CNT files. Only intended for
#' import of 32-bit files.
#'
#' @param file_name Name of .CNT file to be loaded.
#' @importFrom tibble tibble
#' @keywords internal

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
  out <- list(chan_info = chan_df,
              head_info = data_info,
              chan_data = chan_data,
              event_list = ev_list)
}

#' Function for import Brain Vision Analyzer files
#' @param file_name file name of the header file.
#' @keywords internal
import_vhdr <- function(file_name) {

  .data <- read_vhdr(file_name)
  epochs <- data.frame(epoch = 1,
                       recording = basename(tools::file_path_sans_ext(file_name)))
  .data$epochs <- epochs
  .data
}

#' @importFrom ini read.ini
#' @keywords internal
read_vhdr <- function(file_name) {

  header_info <- ini::read.ini(file_name)

  if (identical(header_info$`Common Infos`$DataFormat, "ASCII")) {
    stop("BVA ASCII not currently supported.")
  }

  if (identical(header_info$`Common Infos`$DataType,
                "FREQUENCYDOMAIN")) {
    stop("Only time domain data is currently supported,
     this data is in the frequency domain.")
  }

  if (identical(header_info$`Common Infos`$DataOrientation,
                "VECTORIZED")) {
    multiplexed <- FALSE
  } else {
    multiplexed <- TRUE
  }

  data_file <- header_info$`Common Infos`$DataFile
  vmrk_file <- header_info$`Common Infos`$MarkerFile
  n_chan <- as.numeric(header_info$`Common Infos`$NumberOfChannels)
  # header gives sampling times in microseconds
  srate <- 1e6 / as.numeric(header_info$`Common Infos`$SamplingInterval)
  bin_format <- header_info$`Binary Infos`$BinaryFormat

  chan_labels <- lapply(header_info$`Channel Infos`,
                        function(x) unlist(strsplit(x = x, split = ",")))
  chan_labels <- sapply(chan_labels, "[[", 1)
  chan_info <- parse_vhdr_chans(chan_labels,
                                header_info$`Coordinates`)

  file_size <- file.size(data_file)
  .data <- read_dat(data_file,
                    file_size = file_size,
                    bin_format = bin_format,
                    multiplexed)

  .data <- matrix(.data, ncol = n_chan, byrow = multiplexed)
  .data <- tibble::as.tibble(.data)
  names(.data) <- chan_labels
  n_points <- nrow(.data)
  timings <- tibble::tibble(sample = 1:n_points,
                            time = (sample - 1) / srate)

  .markers <- read_vmrk(vmrk_file)
  .markers <-
    dplyr::mutate(.markers,
                  event_onset = as.numeric(event_onset),
                  length = as.numeric(length),
                  .electrode = ifelse(.electrode == 0,
                                      NA,
                                      as.character(chan_info$electrode[.electrode])),
                  event_time = (event_onset - 1) / srate)
  .data <- eeg_data(data = .data,
                    srate = srate,
                    events = .markers,
                    chan_info = chan_info,
                    timings = timings,
                    reference = NULL,
                    continuous = TRUE
                    )
  .data
}

#' Read the raw data for a BVA file.
#'
#' @param file_name Filename of the .dat file
#' @param file_size Size of the .dat file
#' @param bin_format Storage format.
#' @param multiplexed Logical; whether the data is VECTORIZED (FALSE) or MULTIPLEXED (TRUE)
#' @keywords internal

read_dat <- function(file_name,
                     file_size,
                     bin_format,
                     multiplexed) {

  if (identical(bin_format, "IEEE_FLOAT_32")) {
    nbytes <- 4
    signed <- TRUE
  } else if (identical(bin_format, "UINT_16")) {
    nbytes <- 2
    signed <- FALSE
  } else {
    nbytes <- 2
    signed <- TRUE
  }

  raw_data <- readBin(file_name,
                      what = "double",
                      n = file_size,
                      size = nbytes,
                      signed = signed)
}

#' Read BVA markers
#'
#' Import and parse BVA marker files.
#'
#' @param file_name File name of the .vmrk markers file.
#' @keywords internal
read_vmrk <- function(file_name) {
  vmrks <- ini::read.ini(file_name)
  marker_id <- names(vmrks$`Marker Infos`)

  markers <- lapply(vmrks$`Marker Infos`, function(x) unlist(strsplit(x, ",")))

  # to do - fix to check for "New segment" instead - it's the extra date field that causes issues
  if (length(markers[[1]]) == 6) {
    date <- markers[[1]][[6]]
    markers[[1]] <- markers[[1]][1:5]
  } else {
    date <- NA
  }

  markers <- tibble::as.tibble(do.call(rbind, markers))
  names(markers) <- c("BVA_type",
                      "event_type",
                      "event_onset",
                      "length",
                      ".electrode")
  markers <- dplyr::mutate(markers,
                           event_onset = as.numeric(event_onset),
                           length = as.numeric(length),
                           .electrode = as.numeric(.electrode))

  markers
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
    fdt_file <- paste0(tools::file_path_sans_ext(file_name),
                       ".fdt")
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

  signals <- data.frame(cbind(t(signals),
                              times))
  srate <- c(temp_dat$EEG[[which(var_names == "srate")]])
  names(signals) <- c(unique(chan_info$electrode), "time")
  signals <- dplyr::group_by(signals, time)
  signals <- dplyr::mutate(signals, epoch = 1:dplyr::n())
  signals <- dplyr::ungroup(signals)

  event_info <- temp_dat$EEG[[which(var_names == "event")]]

  event_table <- tibble::as_tibble(t(matrix(as.integer(event_info),
                               nrow = dim(event_info)[1],
                               ncol = dim(event_info)[3])))

  names(event_table) <- unlist(dimnames(event_info)[1])
  event_table$event_time <- (event_table$latency - 1) / srate
  event_table <- event_table[c("latency", "event_time", "type", "epoch")]
  event_table <- dplyr::rename(event_table,
                               event_type = "type",
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
    event_table$time <- timings[which(timings$sample %in% event_table$event_onset,
                                      arr.ind = TRUE), ]$time
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
#' @keywords internal
parse_chaninfo <- function(chan_info) {
  chan_info <- tibble::as_tibble(t(as.data.frame(chan_info)))
  chan_info <- tibble::as_tibble(data.frame(lapply(chan_info,
                                                   unlist),
                                            stringsAsFactors = FALSE))
  chan_info <- chan_info[, sort(names(chan_info))]
  expected <- c("labels", "radius",
                "ref", "sph.phi",
                "sph.radius", "sph.theta",
                "theta", "type",
                "urchan", "X", "Y", "Z")
  if (!all(names(chan_info) == expected)) {
    warning("EEGLAB chan info has unexpected format, taking electrode names only.")
    out <- data.frame(chan_info["labels"])
    names(out) <- "electrode"
  }
  names(chan_info) <- c("electrode",
                        "radius",
                        "ref",
                        "sph_phi",
                        "sph_radius",
                        "sph_theta",
                        "angle",
                        "type",
                        "urchan",
                        "cart_x",
                        "cart_y",
                        "cart_z"
                        )
  chan_info <- chan_info[c("electrode",
                           "cart_x",
                           "cart_y",
                           "cart_z")]
  # EEGLAB co-ordinates are rotated 90 degrees compared to our coordinates,
  # and left-right flipped
  #chan_info <- rotate_sph(chan_info, -90)
  chan_info$cart_y <- -chan_info$cart_y
  names(chan_info) <- names(chan_info)[c(1, 3, 2, 4)]
  chan_info <- chan_info[, c(1, 3, 2, 4)]
  #chan_info <- flip_x(chan_info)
  sph_coords <- cart_to_spherical(chan_info[, c("cart_x", "cart_y", "cart_z")])
  xy <- project_elecs(sph_coords)
  chan_info <- dplyr::bind_cols(electrode = as.character(chan_info$electrode),
                                sph_coords,
                                chan_info[, 2:4],
                                xy)
  chan_info
  #tibble::as.tibble(chan_info)
  }

#' Parse BVA channel info
#'
#' @param chan_labels The `Channel Infos` section from a BVA vhdr file
#' @param chan_info The `Coordinates` section from a BVA vhdr file
#' @keywords internal
parse_vhdr_chans <- function(chan_labels,
                             chan_info) {
  init_chans <- data.frame(electrode = chan_labels)
  coords <- lapply(chan_info,
                   function(x) as.numeric(unlist(strsplit(x, split = ","))))

  new_coords <- data.frame(do.call(rbind, coords))
  names(new_coords) <- c("radius", "theta", "phi")
  new_coords <- cbind(init_chans, new_coords)

  new_coords[new_coords$radius == 0, 2:4] <- NA

  chan_info <- bva_elecs(new_coords)
  tibble::as.tibble(chan_info)
}

#' Convert BVA spherical locations
#'
#' Reads the BrainVoyager spherical electrode locations and converts them to
#' Cartesian 3D and 2D projections.
#'
#' @param chan_info A data.frame containing electrodes
#' @param circumference Head circumference in millimetres
#' @keywords internal
bva_elecs <- function(chan_info, circumference = 85) {
  chan_info <-
    dplyr::mutate(chan_info,
                  cart_x = sin(theta * pi / 180) * cos(phi * pi / 180),
                  cart_y = sin(theta * pi / 180) * sin(phi * pi / 180),
                  cart_z = cos(theta * pi / 180),
                  x = theta * cos(phi * pi / 180),
                  y = theta * sin(phi * pi /180))
  chan_info[, c("cart_x", "cart_y", "cart_z")] <-
    chan_info[, c("cart_x", "cart_y", "cart_z")] * circumference
  chan_info
}


# Export continuous data in Brain Vision Analyzer format
#' @noRd
export_bva <- function(.data,
                       filename,
                       orientation) {
  UseMethod("export_bva", .data)
}

export_bva.default <- function(.data,
                               filename,
                               orientation) {
  stop("export_bva() can currently only export continuous eeg_data objects.")
}

export_bva.eeg_epochs <- function(.data,
                                  filename,
                                  orientation = "VECTORIZED") {

  if (!(orientation %in% c("VECTORIZED", "MULTIPLEXED"))) {
    stop("Orientation must be VECTORIZED or MULTIPLEXED.")
  }

  filename <- tools::file_path_sans_ext(filename)

  write_vhdr(.data, filename, orientation)
  write_dat(.data, filename, orientation)
  write_vmrk(.data, filename)
}


write_vhdr <- function(.data,
                       filename,
                       orientation) {
  vhdr_file <- paste0(filename, ".vhdr")
  con <- file(vhdr_file, open = "w", encoding = "UTF-8")

  on.exit(close(con))
  new_header <- list()
  new_header[["Common Infos"]] <-
    list(Codepage = "UTF-8",
         DataFile = paste0(filename, ".dat"),
         MarkerFile = paste0(filename, ".vmrk"),
         DataFormat = "BINARY",
         DataOrientation=orientation,
         DataType = "TIMEDOMAIN",
         NumberOfChannels = ncol(.data$signals),
         DataPoints = nrow(.data$signals) * ncol(.data$signals),
         SamplingInterval = 1e6 / .data$srate)
  #new_header[["User Infos"]] <- ""
  new_header[["Binary Infos"]] <- list(BinaryFormat = "IEEE_FLOAT_32")
  new_header[["Channel Infos"]] <- as.list(paste(names(.data$signals),
                                                 .data$reference$ref_chans,
                                                 "",
                                                 paste0(intToUtf8(0x03BC), "V"),
                                                 sep = ","))
  names(new_header[["Channel Infos"]]) <- paste0("Ch", 1:ncol(.data$signals))
  new_header[["Coordinates"]] <- as.list(paste(.data$chan_info$sph_radius,
                                       .data$chan_info$sph_theta,
                                       .data$chan_info$sph_phi,
                                       sep = ","))
  names(new_header[["Coordinates"]]) <- paste0("Ch", 1:ncol(.data$signals))
  writeLines("Brain Vision Data Exchange Header File Version 2.0", con)
  writeLines(paste0("; Created using eegUtils http://craddm.github.io/eegUtils/"),
             con)
  writeLines("", con)

  for (i in names(new_header)) {
    writeLines(paste0("[", i,"]"), con)
    entries <-
      lapply(seq_along(new_header[[i]]),
             function(x) paste0(names(new_header[[i]][x]),
                                "=",
                                new_header[[i]][x]))
    lapply(entries,
           writeLines,
           con)
    writeLines("", con)
  }
  writeLines("", con)
}

write_dat <- function(.data,
                      filename,
                      orientation) {
  vdat_file <- paste0(filename, ".dat")
  con <- file(vdat_file, open = "wb")
  on.exit(close(con))
  if (identical(orientation, "VECTORIZED")) {
    # convert to matrix and vector before writiing
    writeBin(as.vector(as.matrix(.data$signals)), con, size = 4)
  } else if (identical(orientation, "MULTIPLEXED")) {
    # transpose matrix before vectorizing
    writeBin(as.vector(t(as.matrix(.data$signals))), con, size = 4)
  }
}

write_vmrk <- function(.data,
                       filename) {
  vmrk_file <- paste0(filename, ".vmrk")
  con <- file(vmrk_file, open = "w")
  on.exit(close(con))
  writeLines("Brain Vision Data Exchange Header File Version 2.0", con)
  writeLines(paste0("; Created using eegUtils http://craddm.github.io/eegUtils/"),
             con)
  writeLines("", con)

  new_markers <- list("Common Infos" = list(Codepage = "UTF-8",
                                            DataFile = paste0(filename, ".dat")),
                      "Marker Infos" = as.list(paste("Stimulus",
                                      .data$events$event_type,
                                      .data$events$event_onset,
                                      1,
                                      0,
                                      sep = ",")))
  names(new_markers[["Marker Infos"]]) <-
    paste0("Mk", 1:nrow(.data$events))

  for (i in names(new_markers)) {
    writeLines(paste0("[", i,"]"), con)
    entries <-
      lapply(seq_along(new_markers[[i]]),
             function(x) paste0(names(new_markers[[i]][x]),
                                "=",
                                new_markers[[i]][x]))
    lapply(entries,
           writeLines,
           con)

    writeLines("", con)
  }
}
