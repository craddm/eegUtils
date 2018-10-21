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
#' @return object of class \code{eeg_data}, re-referenced as requested.
#' @describeIn reref_eeg Rereference objects of class \code{eeg_data}
#' @export

reref_eeg.eeg_data <- function(data,
                               ref_chans = "average",
                               exclude = NULL,
                               robust = FALSE,
                               ...) {

  n_chans <- ncol(data$signals)

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

  # Get excluded channel names and/or convert to numbers if necessary
  if (!is.null(exclude)) {
    if (is.numeric(exclude)) {
      excluded <- names(data$signals)[exclude]
    } else {
        excluded <- exclude
        exclude <- which(names(data$signals) %in% exclude)
    }
  } else {
    excluded <- NULL
  }

  # Calculate new reference data
  if (ref_chans == "average") {
    if (is.null(exclude)) {
      if (robust) {
        ref_data <- matrixStats::rowMedians(as.matrix(data$signals))
      } else {
        ref_data <- rowMeans(data$signals)
      }
    } else {
      if (robust) {
        ref_data <- matrixStats::rowMedians(as.matrix(data$signals[-exclude]))
      } else {
        ref_data <- rowMeans(data$signals[-exclude])
      }
    }
    #remove reference from data
    data$signals <- data$signals - ref_data
  } else {
    if (any(all(ref_chans %in% colnames(data$signals)) | is.numeric(ref_chans))) {
      if (length(ref_chans) > 1) {
        ref_data <- rowMeans(data$signals[, ref_chans])
      } else {
        ref_data <- unlist(data$signals[, ref_chans])
      }
      data$signals <- data$signals - ref_data
    } else {
      stop("Electrode(s) not found.")
    }
  }
  data$signals <- as.data.frame(data$signals)
  if (ref_chans == "average") {
    data$reference <- list(ref_chans = ref_chans,
                           excluded = excluded)
    } else {
      data$reference <- list(ref_chans = ref_chans,
                             excluded = excluded)
      data <- select_elecs(data, ref_chans, keep = FALSE)
      }
  return(data)
}

#' Baseline correction
#'
#' Used to remove the mean of a specified time period from the data. Currently
#' only performs subtractive baseline. With a data frame, searches for
#' "electrode" and "epoch" columns, and groups on these when found. An electrode
#' column is always required; an epoch column is not.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Data to be baseline corrected.
#' @param ... other parameters to be passed to functions
#' @export

rm_baseline <- function(data, ...) {
  UseMethod("rm_baseline", data)
}

#' @param time_lim Numeric character vector (e.g. time_lim <- c(-.1, 0)). If
#'   none given, defaults to mean of the whole of each epoch if the data is epoched, or the
#'   channel mean if the data is continuous.
#' @describeIn rm_baseline remove baseline from continuous \code{eeg_data}
#' @export

rm_baseline.eeg_data <- function(data, time_lim = NULL, ...) {

  if (is.null(time_lim)) {
    baseline_dat <- colMeans(data$signals)
  } else {
    base_times <- select_times(data, time_lim = time_lim)
    baseline_dat <- colMeans(base_times$signals)
  }
  data$signals <- sweep(data$signals, 2, baseline_dat, '-')
  return(data)
}

#' @describeIn rm_baseline Remove baseline from eeg_epochs
#' @export

rm_baseline.eeg_epochs <- function(data,
                                   time_lim = NULL,
                                   ...) {

  n_epochs <- length(unique(data$timings$epoch))
  n_times <- length(unique(data$timings$time))
  n_chans <- ncol(data$signals)
  elecs <- names(data$signals)

  if (is.null(time_lim)) {
    data$signals <- as.matrix(data$signals)
    dim(data$signals) <- c(n_times, n_epochs, n_chans)
    baseline_dat <- colMeans(data$signals)
    data$signals <- vapply(1:n_times,
                           function(x) data$signals[x, , ] - baseline_dat,
                           array(0, dim = c(n_epochs, n_chans)))
    data$signals <- aperm(data$signals, c(3, 1, 2))
    data$signals <- array(data$signals, dim = c(n_epochs * n_times, n_chans))
    data$signals <- as.data.frame(data$signals)
    names(data$signals) <- elecs
  } else {
    base_times <- select_times(data,
                               time_lim = time_lim)
    base_times$signals <- split(base_times$signals,
                                base_times$timings$epoch)
    base_baselines <- lapply(base_times$signals,
                             Matrix::colMeans)
    data$signals <- split(data$signals,
                          data$timings$epoch)
    data$signals <- lapply(seq_along(data$signals),
                           function(i) sweep(data$signals[[i]],
                                             2,
                                             base_baselines[[i]],
                                             "-"))
    data$signals <- do.call(rbind, data$signals)
  }
  data
}

#' @describeIn rm_baseline Legacy method for data.frames
#' @export
rm_baseline.data.frame <- function(data,
                                   time_lim = NULL,
                                   ...) {

  if (!("time" %in% colnames(data))) {
    stop("Time dimension is required.")
  }

  if (length(time_lim) == 1) {
    stop("time_lim should specify the full time range.")
  }

  # if the data is epoched, group by electrode and epoch; otherwise, just by
  # electrode.

  if ("epoch" %in% colnames(data)) {
    data <- dplyr::group_by(data,
                            electrode,
                            epoch,
                            add = TRUE)
  } else{
    data <- dplyr::group_by(data,
                            electrode,
                            add = TRUE)
  }

  if (is.null(time_lim)) {
    # if no time_lim provided, just delete mean of all time points
    data <- dplyr::mutate(data,
                          amplitude = amplitude - mean(amplitude))
  } else {

    data_sel <- dplyr::filter(data,
                              time >= time_lim[1],
                              time <= time_lim[2])
    baseline <- dplyr::summarise(data_sel,
                                 bl = mean(amplitude))
    # This is relatively memory intensive - not so bad now but would prefer
    # another way. Could get extremely painful with time-frequency data.
    data <- dplyr::left_join(data,
                             baseline)
    data <- dplyr::mutate(data,
                          amplitude = amplitude - bl)
    data <- dplyr::select(data,
                          -bl)
  }
  data <- ungroup(data)
  data
}

#' @param type Type of baseline correction to apply. Options are ("divide",
#'   "ratio", "absolute", "db")
#' @describeIn rm_baseline Method for \code{eeg_tfr} objects
rm_baseline.eeg_tfr <- function(data,
                                time_lim = NULL,
                                type = "divide",
                                ...) {
  bline <- select_times(data, time_lim)
  bline <- apply(bline$signals,
                 c(2, 3),
                 mean,
                 na.rm = TRUE)

  do_corrs <- function(data,
                       type,
                       bline) {
    switch(type,
           "divide" = ((data - bline) / bline) * 100,
           "absolute" = data - bline,
           "db" = 10 * log10(data / bline),
           "ratio" = data / bline,
           )
  }
  zz <- apply(data$signals,
              1,
              do_corrs,
              type = type,
              bline = bline)

              #function(x) x / bline)
  orig_dims <- dim(data$signals)
  dim(zz) <- c(orig_dims[2],
               orig_dims[3],
               orig_dims[1])
  zz <- aperm(zz,
              c(3, 1, 2))
  data$signals[,,] <- zz
  data$freq_info$baseline <- type
  data
}

#' Create epochs from EEG data
#'
#' Creates epochs around specified event triggers. Requires data of class
#' \code{eeg_data}. Where multiple events are specified, epochs will be created
#' around each event.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Continuous data to be epoched.
#' @param ... Parameters passed to functions
#' @export

epoch_data <- function(data, ...) {
  UseMethod("epoch_data", data)
}

#' Create epochs from EEG data
#'
#' @param data Continuous data to be epoched.
#' @param ... Parameters passed to functions
#' @export

epoch_data.default <- function(data, ...) {
  stop("Requires object of class eeg_data.")
}

#' @param events Character vector of events to epoch around.
#' @param time_lim Time in seconds to form epoch around the events. Defaults to
#'   one second either side.
#' @importFrom dplyr left_join
#' @importFrom purrr map map_df
#'
#' @return Returns an epoched object of class \code{eeg_data}
#'
#' @describeIn epoch_data
#'
#' @export
#'

epoch_data.eeg_data <- function(data,
                                events,
                                time_lim = c(-1, 1),
                                ...) {

  if (!any(events %in% unique(data$events$event_type))) {
    stop("No events found - check event codes.")
  }

  if (!all(events %in% unique(data$events$event_type))) {
    warning("Some events not found - check event codes.")
  }

  # If the data has been downsampled, sample spacing will be greater than 1.
  # Subsequent steps need to account for this when selecting based on sample number.
  # Multiplying srate by spacing accomplishes this.

  samp_diff <- unique(diff(data$timings$sample))
  if (length(samp_diff) > 1) {
    stop("Sample spacing is uneven, cannot downsample.")
  } else {
    srate <- data$srate * samp_diff
  }

  # create a vector that counts back and ahead of the timelock event in samples
  # i.e. if srate is 1000, a vector from -100 to 0 to 100 would be -.1 s before
  # to .1 s after event onset

  samps <- seq(round(time_lim[[1]] * srate),
               round(time_lim[[2]] * (srate - 1)),
               by = samp_diff)

  event_table <- data$events

  # go through all events and find which match those specified as the events to
  # epoch around
  epoch_zero <-
    sort(unlist(purrr::map(events,
                           ~ event_table[which(event_table$event_type == .), ]$event_onset)))

  # for each epoching event, create a vector of samples before and after that
  # are within the time limits
  epoched_data <- purrr::map(epoch_zero,
                             ~ . + samps)

  epoched_data <- purrr::map_df(epoched_data,
                                ~ tibble::tibble(sample = .,
                                                 time = samps / srate),
                                .id = "epoch")

  epoched_data$epoch <- as.numeric(epoched_data$epoch)

  # create new event_table
  event_table <- dplyr::inner_join(event_table,
                                   epoched_data,
                                   by = c("event_onset" = "sample"))

  epoched_data <- dplyr::left_join(epoched_data,
                                   cbind(data$signals,
                                         data$timings),
                                   by = c("sample" = "sample"))

  # Check for any epochs that contain NAs
  epoched_data <- split(epoched_data, epoched_data$epoch)
  na_epochs <- vapply(epoched_data,
                      function(x) any(is.na(x)),
                      FUN.VALUE = logical(1))
  epoched_data <- epoched_data[!na_epochs]
  epoched_data <- do.call("rbind", epoched_data)

  event_table <- event_table[event_table$epoch %in% names(na_epochs[!na_epochs]), ]

  if (any(na_epochs)) {
    message(paste(sum(na_epochs),
                  "epoch(s) with NAs removed. Epoch boundaries would lie outside data."))
  }

  epoched_data$time.y <- NULL
  names(epoched_data)[[3]] <- "time"
  data$signals <- epoched_data[, -1:-3] #check which columns this is...
  data$timings <- epoched_data[, 1:3]
  data$continuous <- FALSE
  data$events <- event_table
  class(data) <- c("eeg_epochs", "eeg_data")
  return(data)
}

#' Create epochs
#'
#' @param data Continuous data to be epoched
#' @param ... Parameters passed to methods
#' @export
#'

epoch_data.eeg_epochs <- function(data, ...) {
  stop("Data is already epoched, cannot currently re-epoch.")
}

#' Downsampling EEG data
#'
#' Performs low-pass anti-aliasing filtering and downsamples EEG data by a
#' specified factor. This is a wrapper for \code{decimate} from the
#' \code{signal} package. Note that this will also adjust the event table,
#' moving events to the nearest time remaining after downsampling
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

eeg_downsample.eeg_epochs <- function(data, q, ...) {

  q <- as.integer(q)

  if (q < 2) {
    stop("q must be 2 or more.")
  } else if ((data$srate / q) %% 1 > 0){
    stop("srate / q must give a round number.")
  }

  message(paste0("Downsampling from ", data$srate, "Hz to ",
                 data$srate / q, "Hz."))

  epo_length <- length(unique(data$timings$time)) %% q

  #ceiling(epo_length / q)

  if (epo_length > 0) {
    message("Dropping ", epo_length, " time points to make n of samples a multiple of q.")
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

