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
#' reref_eeg(demo_epochs, exclude = "B5")
#' # Reference data using the median of the reference channels rather than the mean
#' reref_eeg(demo_epochs, robust = TRUE)
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to re-reference. Primarily meant for use with data of class
#'   \code{eeg_data}.
#' @param ... Further parameters to be passed to reref_eeg
#' @export

reref_eeg <- function(data, ...) {
  UseMethod("reref_eeg", data)
}

#' @export
#' @describeIn reref_eeg Default method
reref_eeg.default <- function(data, ...) {
  stop(paste("reref_eeg does not know how to handle data of class", class(data)))
}

#' @param ref_chans Channels to reference data to. Defaults to "average" i.e.
#'   average of all electrodes in data. Character vector of channel names or numbers.
#' @param exclude Electrodes to exclude from average reference calculation.
#' @param robust Use median instead of mean; only used for average reference.
#' @importFrom matrixStats rowMedians
#' @import data.table
#' @return object of class \code{eeg_data}, re-referenced as requested.
#' @describeIn reref_eeg Rereference objects of class \code{eeg_data}
#' @export

reref_eeg.eeg_data <- function(data,
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
      data <- select_elecs(data, ref_chans, keep = FALSE)
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

eeg_downsample.eeg_data <- function(data, q, ...) {

  q <- as.integer(q)

  if (q < 2) {
    stop("q must be 2 or more.")
  } else if ((data$srate / q) %% 1 > 0){
    stop("srate / q must give a round number.")
  }

  message(paste0("Downsampling from ", data$srate, "Hz to ",
                 data$srate / q, "Hz."))

  data_length <- length(unique(data$timings$time)) %% q

  #ceiling(epo_length / q)

  if (data_length > 0) {
    message("Dropping ",
            data_length,
            " time points to make n of samples a multiple of q.")
    new_times <- utils::head(unique(data$timings$time),
                      -data_length)
    data <- select_times(data,
                         time_lim = c(min(new_times),
                                      max(new_times)))
  }

  # step through each column and decimate each channel
  data$signals <- purrr::map_df(as.list(data$signals),
                                ~signal::decimate(., q))

  # select every qth timing point, and divide srate by q
  data$srate <- data$srate / q
  data$timings <- data$timings[seq(1, length(data$timings[[1]]), by = q), ]

  # The event table also needs to be adjusted. Note that this inevitably jitters
  # event timings by up to q/2 sampling points.
  nearest_samps <- findInterval(data$events$event_onset,
                                data$timings$sample)
  data$events$event_onset <- data$timings$sample[nearest_samps]
  data$events$event_time <- data$timings$time[nearest_samps]
  data
}

#' @describeIn eeg_downsample Downsample eeg_epochs objects
#' @export

eeg_downsample.eeg_epochs <- function(data,
                                      q,
                                      ...) {

  q <- as.integer(q)

  if (q < 2) {
    stop("q must be 2 or more.")
  } else if ((data$srate / q) %% 1 > 0){
    stop("srate / q must give a round number.")
  }

  message(paste0("Downsampling from ",
                 data$srate, "Hz to ",
                 data$srate / q, "Hz."))

  epo_length <- length(unique(data$timings$time)) %% q

  if (epo_length > 0) {
    message("Dropping ",
            epo_length,
            " time points to make n of samples a multiple of q.")
    new_times <- utils::head(unique(data$timings$time),
                             -epo_length)
    data <- select_times(data,
                         time_lim = c(min(new_times),
                                      max(new_times)))
  }

  data$signals <- split(data$signals, data$timings$epoch)

  new_times <- data$timings$time
  new_length <- nrow(data$signals[[1]]) #- epo_length
  data$signals <- lapply(data$signals,
                         `[`,
                         1:new_length,
                         )
  data$signals <- purrr::map_df(data$signals,
                             ~purrr::map_df(as.list(.),
                                            ~signal::decimate(., q)))
  # step through each column and decimate each channel
  # data$signals <- purrr::map_df(as.list(data$signals),
  #                               ~signal::decimate(., q))

  # select every qth timing point, and divide srate by q
  data$srate <- data$srate / q
  data$timings <- data$timings[seq(1,
                                   length(data$timings[[1]]),
                                   by = q), ]

  # The event table also needs to be adjusted. Note that this inevitably jitters
  # event timings by up to q/2 sampling points.
  data_samps <- sort(unique(data$timings$sample))
  #samp_times <-
  nearest_samps <- findInterval(data$events$event_onset,
                                data_samps)
  data$events$event_onset <- data_samps[nearest_samps]
  data$events$event_time <- 1 / (data$srate * q) * data$events$event_onset
  data
}

#' Combine EEG objects
#'
#' Combine multiple \code{eeg_epochs} or \code{eeg_data} objects into a single
#' object. Note that this does not currently perform any sort of checking for
#' internal consistency or duplication. It simply combines the objects in the
#' order they are passed.
#'
#' @param data An \code{eeg_epochs} object
#' @param ... additional \code{eeg_epochs} objects
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export
#'

eeg_combine <- function(data, ...) {
  UseMethod("eeg_combine", data)
}

#'@export
eeg_combine.default <- function(data, ...) {
  stop(
    "Don't know how to combine objects of class",
    class(data)[[1]],
    call. = FALSE
  )
}

#' @describeIn eeg_combine Method for combining \code{eeg_data} objects.
#' @export

eeg_combine.eeg_data <- function(data, ...){

  args <- list(...)
  if (length(args) == 0) {
    stop("Nothing to combine.")
  }
  if (all(sapply(args, is.eeg_data))) {
    if (any(sapply(args, is.eeg_epochs))) {
      stop("All objects must be unepoched eeg_data objects.")
    } else {
      data$signals <- dplyr::bind_rows(data$signals,
                                       purrr::map_df(args, ~.$signals))
      data$events <- dplyr::bind_rows(data$events,
                                      purrr::map_df(args, ~.$events))
      data$timings <- dplyr::bind_rows(data$timings,
                                       purrr::map_df(args, ~.$timings))
    }
  }
  data
}

#' @describeIn eeg_combine Method for combining \code{eeg_epochs} objects
#' @export

eeg_combine.eeg_epochs <- function(data, ...) {

  args <- list(...)
  if (length(args) == 0) {
    stop("Nothing to combine.")
  }
  if (all(sapply(args, is.eeg_epochs))) {
    data$signals <- dplyr::bind_rows(data$signals,
                                     purrr::map_df(args, ~.$signals))
    data$events <- dplyr::bind_rows(data$events,
                                    purrr::map_df(args, ~.$events))
    data$timings <- dplyr::bind_rows(data$timings,
                        purrr::map_df(args, ~.$timings))
  } else {
    stop("All inputs must be eeg_epochs objects.")
  }
  #fix epoch numbering for combined objects
  data <- check_timings(data)
  data
}


#' Check consistency of event and timing tables
#'
#' @param data \code{eeg_epochs} object
#' @keywords internal

check_timings <- function(data) {

  n_rows <- nrow(data$timings)
  epochs <- unique(data$timings$epoch)
  duplicated()

  # if the epoch numbers are not ascending, fix them...
  while (any(diff(data$timings$epoch) < 0)) {

    # check consistency of the data timings table.
    # timings should be consistently increasing.
    # only works correctly with 2 objects

    #check for any places where epoch numbers decrease instead of increase
    switch_locs <- which(diff(data$timings$epoch) == min(diff(data$timings$epoch)))

    #consider switch this out with an RLE method, which would be much simpler.

    if (length(switch_locs) == 1) {
      switch_epo <- data$timings$epoch[switch_locs]
      switch_sample <- data$timings$sample[switch_locs]
      data$timings$epoch[(switch_locs + 1):n_rows] <-
        data$timings$epoch[(switch_locs + 1):n_rows] + switch_epo
      data$timings$sample[(switch_locs + 1):n_rows] <-
        data$timings$sample[(switch_locs + 1):n_rows] + switch_sample
    } else {
      for (i in 1:(length(switch_locs) - 1)) {
        switch_epo <- data$timings$epoch[switch_locs[i]]
        switch_sample <- data$timings$sample[switch_locs[i]]
        data$timings$epoch[(switch_locs[i] + 1):switch_locs[i + 1]] <-
          data$timings$epoch[(switch_locs[i] + 1):switch_locs[i + 1]] + switch_epo
        data$timings$sample[(switch_locs[i] + 1):switch_locs[i + 1]] <-
          data$timings$sample[(switch_locs[i] + 1):switch_locs[i + 1]] + switch_sample
      }

      switch_epo <- data$timings$epoch[switch_locs[length(switch_locs)]]
      switch_sample <- data$timings$sample[switch_locs[length(switch_locs)]]
      data$timings$epoch[(switch_locs[length(switch_locs)] + 1):n_rows] <-
        data$timings$epoch[(switch_locs[length(switch_locs)] + 1):n_rows] + switch_epo
      data$timings$sample[(switch_locs[length(switch_locs)] + 1):n_rows] <-
        data$timings$sample[(switch_locs[length(switch_locs)] + 1):n_rows] + switch_sample
    }
  }

  if (any(diff(data$events$event_time) < 0)) {

    #check consistency of the data event table
    #to handle cases where there are multiple events per epoch,
    #use RLE to code how many times each epoch number comes up in a row;
    #then replace old epoch numbers with new and reconstruct the vector
    orig_ev <- rle(data$events$epoch)
    orig_ev$values <- unique(data$timings$epoch)
    data$events$epoch <- inverse.rle(orig_ev)
    data$events <- dplyr::left_join(data$events,
                                    data$timings,
                                    by = c("epoch", "time"))
    data$events$event_onset <- data$events$sample
    data$events$sample <- NULL
  }
  data
}

