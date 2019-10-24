#' Referencing
#'
#' Used to reference the EEG data to a specified electrode or electrodes.
#' Defaults to average reference. When specific electrodes are used, they are
#' removed from the data. Meta-data about the referencing scheme is held in the
#' \code{eeg_data} structure.
#'
#' @examples
#' # demo_epochs is average referenced by default
#' demo_epochs
#' # Rereference it but exclude B5 from calculation of the average
#' eeg_reference(demo_epochs, exclude = "B5")
#' # Reference data using the median of the reference channels rather than the mean
#' eeg_reference(demo_epochs, robust = TRUE)
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to re-reference. Primarily meant for use with data of class
#'   \code{eeg_data}.
#' @param ... Further parameters to be passed to eeg_reference
#' @export

eeg_reference <- function(data, ...) {
  UseMethod("eeg_reference", data)
  }

#' @export
#' @describeIn eeg_reference Default method
eeg_reference.default <- function(data, ...) {
  stop("eeg_reference does not know how to handle data of class ",
       paste0(class(data), collapse = " "))
}

#' @param ref_chans Channels to reference data to. Defaults to "average" i.e.
#'   average of all electrodes in data. Character vector of channel names or
#'   numbers.
#' @param exclude Electrodes to exclude from average reference calculation.
#' @param robust Use median instead of mean; only used for average reference.
#' @importFrom matrixStats rowMedians
#' @import data.table
#' @return object of class \code{eeg_data}, re-referenced as requested.
#' @describeIn eeg_reference Rereference objects of class \code{eeg_data}
#' @export

eeg_reference.eeg_data <- function(data,
                                   ref_chans = "average",
                                   exclude = NULL,
                                   robust = FALSE,
                                   ...) {

  # check for existing reference. Add it back in if it exists.
  if (!is.null(data$reference)) {
    if (!identical(data$reference$ref_chans, "average")) {
      data$signals[data$reference$ref_chans] <- 0
    }
  }

  # Convert ref_chan channel numbers into channel names
  if (is.numeric(ref_chans)) {
    ref_chans <- names(data$signals)[ref_chans]
  }

  # If average reference is requested, first get all channel names.
  if (identical(ref_chans, "average")) {
    reference <- names(data$signals)
  }

  # Get excluded channel names and/or convert to numbers if necessary
  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      exclude <- names(data$signals)[exclude]
    } else {
      exclude <- exclude[which(exclude %in% names(data$signals))]
    }
    reference <- reference[!(reference %in% exclude)]
  }

  # Calculate new reference data
  if (ref_chans == "average") {
    if (robust) {
      ref_data <- matrixStats::rowMedians(as.matrix(data$signals[, reference]))
    } else {
      ref_data <- rowMeans(data$signals[, reference])
    }
    #remove reference from data
    data$signals <- data.table(data$signals)
    data$signals <- data$signals[, lapply(.SD, function(x) x - ref_data)]
  } else {
    if (any(all(ref_chans %in% colnames(data$signals)) | is.numeric(ref_chans))) {
      if (length(ref_chans) > 1) {
        ref_data <- rowMeans(data$signals[, ref_chans])
      } else {
        ref_data <- unlist(data$signals[, ref_chans])
      }
      data$signals <- data.table(data$signals)
      data$signals <- data$signals[, lapply(.SD, function(x) x - ref_data)]
    } else {
      stop("Electrode(s) not found.")
    }
  }
  data$signals <- tibble::as_tibble(data$signals)
  if (ref_chans == "average") {
    data$reference <- list(ref_chans = ref_chans,
                           excluded = exclude)
    } else {
      data$reference <- list(ref_chans = ref_chans,
                             excluded = exclude)
      data <- select_elecs(data,
                           ref_chans,
                           keep = FALSE)
  }
  data
}

eeg_reference.eeg_ICA <- function(data,
                                  ...) {
  stop("Cannot rereference ICA decompositions.")
}

#' @describeIn eeg_reference Rereference objects of class \code{eeg_data}
#' @export
eeg_reference.eeg_epochs <- function(data,
                                     ref_chans = "average",
                                     exclude = NULL,
                                     robust = FALSE,
                                     ...) {
  do_referencing(data,
                 ref_chans = "average",
                 exclude = NULL,
                 robust = FALSE,
                 ...)
}

do_referencing <- function(data,
                           ref_chans = "average",
                           exclude = NULL,
                           robust = FALSE) {

  # check for existing reference. Add it back in if it exists.
  if (!is.null(data$reference)) {
    if (!identical(data$reference$ref_chans, "average")) {
      data$signals[data$reference$ref_chans] <- 0
    }
  }

  # Convert ref_chan channel numbers into channel names
  if (is.numeric(ref_chans)) {
    ref_chans <- names(data$signals)[ref_chans]
  }

  # If average reference is requested, first get all channel names.
  if (identical(ref_chans, "average")) {
    reference <- names(data$signals)
  }

  # Get excluded channel names and/or convert to numbers if necessary
  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      exclude <- names(data$signals)[exclude]
    } else {
      exclude <- exclude[which(exclude %in% names(data$signals))]
    }
  reference <- reference[!(reference %in% exclude)]
  }

  # Calculate new reference data
  if (identical(ref_chans, "average")) {
    if (robust) {
      ref_data <- matrixStats::rowMedians(as.matrix(data$signals[, reference]))
    } else {
      ref_data <- rowMeans(data$signals[, reference])
    }
  #remove reference from data
    data$signals <- data.table(data$signals)
    data$signals <- data$signals[, lapply(.SD, function(x) x - ref_data)]
  } else {
    if (any(all(ref_chans %in% colnames(data$signals)) | is.numeric(ref_chans))) {
      if (length(ref_chans) > 1) {
        ref_data <- rowMeans(data$signals[, ref_chans])
      } else {
        ref_data <- unlist(data$signals[, ref_chans])
      }
      data$signals <- data.table(data$signals)
      data$signals <- data$signals[, lapply(.SD, function(x) x - ref_data)]
    } else {
      stop("Electrode(s) not found.")
    }
  }

  data$signals <- tibble::as_tibble(data$signals)
  if (identical(ref_chans, "average")) {
    data$reference <- list(ref_chans = ref_chans,
                           excluded = exclude)
    } else {
      data$reference <- list(ref_chans = ref_chans,
                             excluded = exclude)
      data <- select_elecs(data,
                           ref_chans,
                           keep = FALSE)
    }
  data
}

#' Downsampling EEG data
#'
#' Performs low-pass anti-aliasing filtering and downsamples EEG data by a
#' specified factor. This is a wrapper for \code{decimate} from the
#' \code{signal} package. Note that this will also adjust the event table,
#' moving events to the nearest time remaining after downsampling
#'
#' @examples
#' eeg_downsample(demo_epochs, 2)
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An \code{eeg_data} object to be downsampled
#' @param ... Parameters passed to functions
#' @export

eeg_downsample <- function(data, ...) {
  UseMethod("eeg_downsample", data)
}

#' @export
eeg_downsample.default <- function(data, ...) {
  stop("Only used for eeg_data objects at present.")
}

#' @param q Integer factor to downsample by
#' @importFrom signal decimate
#' @describeIn eeg_downsample Downsample eeg_data objects
#' @export

eeg_downsample.eeg_data <- function(data,
                                    q,
                                    ...) {

  q <- check_q(q,
               data$srate)

  data_length <- length(unique(data$timings$time)) %% q

  if (data_length > 0) {
    data <- drop_points(data, data_length)
  }

  # step through each column and decimate each channel
  data$signals <- purrr::map_df(as.list(data$signals),
                                ~signal::decimate(., q))

  # select every qth timing point, and divide srate by q
  data$srate <- data$srate / q
  data$timings <- data$timings[seq(1, length(data$timings[[1]]), by = q), ]

  # The event table also needs to be adjusted. Note that this inevitably jitters
  # event timings by up to q/2 sampling points.
  events(data) <- downsample_events(data$timings,
                                    data$events,
                                    data$srate,
                                    q)
  data
}

#' @describeIn eeg_downsample Downsample eeg_epochs objects
#' @export
eeg_downsample.eeg_epochs <- function(data,
                                      q,
                                      ...) {

  q <- check_q(q,
               data$srate)

  epo_length <- length(unique(data$timings$time)) %% q

  if (epo_length > 0) {
    data <- drop_points(data, epo_length)
  }

  data$signals <- split(data$signals,
                        data$timings$epoch)

  new_times <- data$timings$time
  new_length <- nrow(data$signals[[1]]) #- epo_length
  data$signals <- lapply(data$signals,
                         `[`,
                         1:new_length,
                         )
  # step through each column and decimate each channel, by epoch
  data$signals <-
    purrr::map_df(data$signals,
                  ~purrr::map_df(as.list(.),
                                 ~signal::decimate(., q)))

  # select every qth timing point, and divide srate by q
  data$srate <- data$srate / q
  data$timings <- data$timings[seq(1,
                                   length(data$timings[[1]]),
                                   by = q), ]

  # The event table also needs to be adjusted. Note that this inevitably jitters
  # event timings by up to q/2 sampling points.
  events(data) <- downsample_events(data$timings,
                                    data$events,
                                    data$srate,
                                    q)
  data
}

#' Drop points before downsampling
#'
#' @param data data to be downsampled
#' @param data_length number of points to drop
#' @keywords internal
drop_points <- function(data, data_length) {

  message("Dropping ",
          data_length,
          " time points to make n of samples a multiple of q.")
  new_times <- utils::head(unique(data$timings$time),
                           -data_length)
  # Use custom filter method instead of select_times.
  data <- dplyr::filter(data,
                        time >= min(new_times),
                        time <= max(new_times))
  data
}

#' Downsample the events table
#'
#' @author Matt Craddock \email{matt@@craddock.com}
#' @param timings the timings from the data
#' @param events the events table to downsample
#' @param srate sampling rate
#' @param q downsampling factor
#' @keywords internal
downsample_events <- function(timings,
                              events,
                              srate,
                              q) {

  data_samps <- sort(unique(timings$sample))
  nearest_samps <- findInterval(events$event_onset,
                                data_samps)
  events$event_onset <- data_samps[nearest_samps]
  events$event_time <- 1 / (srate * q) * (events$event_onset - 1)
  events
}

#' Validate the q factor for downsampling
#'
#' @param q Q factor
#' @param srate Sampling rate
#' @keywords internal
check_q <- function(q,
                    srate) {
  q <- as.integer(q)

  if (q < 2) {
    stop("q must be 2 or more.")
  } else if ((srate / q) %% 1 > 0){
      stop("srate / q must give a round number.")
  }

  message(paste0("Downsampling from ",
                 srate, "Hz to ",
                 srate / q, "Hz."))
  q
}



