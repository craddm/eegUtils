#' Function for importing Brain Vision Analyzer files
#' @param file_name file name of the header file.
#' @keywords internal
import_vhdr <- function(file_name,
                        participant_id,
                        recording) {

  .data <- read_vhdr(file_name)
  epochs <- tibble::new_tibble(list(epoch = 1,
                                    participant_id = participant_id,
                                    recording = recording),
                               nrow = 1,
                               class = "epoch_info")
  .data$epochs <- epochs
  .data
}

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

  data_file <- file.path(dirname(file_name),
                         header_info$`Common Infos`$DataFile)
  vmrk_file <- file.path(dirname(file_name),
                         header_info$`Common Infos`$MarkerFile)
  n_chan <- as.numeric(header_info$`Common Infos`$NumberOfChannels)
  # header gives sampling times in microseconds
  srate <- 1e6 / as.numeric(header_info$`Common Infos`$SamplingInterval)
  bin_format <- header_info$`Binary Infos`$BinaryFormat

  chan_labels <- lapply(header_info$`Channel Infos`,
                        function(x) unlist(strsplit(x = x, split = ",")))
  chan_names <- sapply(chan_labels, "[[", 1)
  #EEGLAB exported BVA files are missing the channel resolution, which occupies
  #the 3rd position - check length of chan_labels elements first to cope with
  #that. Probably always 1 microvolt anyway...
  if (length(chan_labels[[1]]) == 3) {
    chan_scale <- as.numeric(sapply(chan_labels, "[[", 3))
  } else {
    chan_scale <- rep(1, length(chan_labels))
  }

  chan_info <- parse_vhdr_chans(chan_names,
                                header_info$`Coordinates`)

  file_size <- file.size(data_file)
  .data <- read_dat(data_file,
                    file_size = file_size,
                    bin_format = bin_format,
                    multiplexed)

  .data <- matrix(.data,
                  ncol = n_chan,
                  byrow = multiplexed)

  if (any(!is.na(chan_scale))) {
    .data <- sweep(.data,
                   2,
                   chan_scale,
                   "*")
  }
  colnames(.data) <- chan_names
  .data <- tibble::as_tibble(.data)

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
                    reference = NULL)
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
    file_type <- "double"
  } else if (identical(bin_format, "UINT_16")) {
    nbytes <- 2
    signed <- FALSE
    file_type <- "int"
  } else {
    nbytes <- 2
    signed <- TRUE
    file_type <- "int"
  }

  raw_data <- readBin(file_name,
                      what = file_type,
                      n = file_size,
                      size = nbytes,
                      signed = signed)
  raw_data
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

  markers <- lapply(vmrks$`Marker Infos`,
                    function(x) unlist(strsplit(x, ",")))

  # to do - fix to check for "New segment" instead - it's the extra date field that causes issues
  if (length(markers[[1]]) == 6) {
    date <- markers[[1]][[6]]
    markers <- lapply(markers, `[`, 1:5)
  } else {
    date <- NA
  }

  markers <- do.call(rbind, markers)

  colnames(markers) <- c("BVA_type",
                         "event_type",
                         "event_onset",
                         "length",
                         ".electrode")

  markers <- tibble::as_tibble(markers)
  markers <-
    dplyr::mutate(markers,
                  event_onset = as.numeric(event_onset),
                  length = as.numeric(length),
                  .electrode = as.numeric(.electrode))
  markers
}

#' Parse BVA channel info
#'
#' @param chan_labels The `Channel Infos` section from a BVA vhdr file
#' @param chan_info The `Coordinates` section from a BVA vhdr file
#' @keywords internal
parse_vhdr_chans <- function(chan_labels,
                             chan_info,
                             verbose = TRUE) {

  init_chans <- data.frame(electrode = chan_labels)
  if (is.null(chan_info)) {
    init_chans$radius <- NA
    init_chans$theta <- NA
    init_chans$phi <- NA
    if (verbose) message("No channel locations found.")
    return(validate_channels(init_chans))
  } else {
    coords <- lapply(chan_info,
                     function(x) as.numeric(unlist(strsplit(x, split = ","))))
    new_coords <- data.frame(do.call(rbind,
                                     coords))
    names(new_coords) <- c("radius", "theta", "phi")
    new_coords <- cbind(init_chans, new_coords)
    new_coords[new_coords$radius %in% c(0, NaN), 2:4] <- NA

    chan_info <- bva_elecs(new_coords)
    tibble::as_tibble(chan_info)
  }
}

#' Convert BVA spherical locations
#'
#' Reads the BrainVoyager spherical electrode locations and converts them to
#' Cartesian 3D and 2D projections.
#'
#' @param chan_info A data.frame containing electrodes
#' @param radius Head radius in millimetres
#' @keywords internal
bva_elecs <- function(chan_info, radius = 85) {
  chan_info <-
    dplyr::mutate(chan_info,
                  cart_x = sin(deg2rad(theta)) * cos(deg2rad(phi)),
                  cart_y = sin(deg2rad(theta)) * sin(deg2rad(phi)),
                  cart_z = cos(deg2rad(theta)),
                  x = theta * cos(deg2rad(phi)),
                  y = theta * sin(deg2rad(phi)))
  chan_info[, c("cart_x", "cart_y", "cart_z")] <-
    chan_info[, c("cart_x", "cart_y", "cart_z")] * radius
  chan_info
}
