#' Referencing
#'
#' Used to reference the EEG data to specified electrode or electrodes. Defaults
#' to average reference. Note that this stores the reference channel in the
#' \code{eeg_data} structure; if one already exists, it is added back to the
#' data before further referencing takes place.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to re-reference. Primarily meant for use with data of class
#'   \code{eeg_data}. Long format is expected if a data-frame is submitted.
#' @param ref_chans Channels to reference data to. Defaults to "average" i.e.
#'   average of all electrodes in data. Character vector of channel names or numbers.
#' @param exclude Electrodes to exclude from average reference calculation.
#' @param robust Use median instead of mean; only applied for average reference.
#' @importFrom tidyr spread
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom matrixStats rowMedians
#' @return object of class \code{eeg_data}, re-referenced as requested.
#' @export
#'

reref_eeg <- function(data, ref_chans = "average", exclude = NULL,
                      robust = FALSE) {

  if (is.eeg_data(data)) {
    tmp_data <- data$signals
    n_chans <- dim(tmp_data)[[2]]
  } else {
    tmp_data <- tidyr::spread(data, electrode, amplitude)
    tmp_data <- dplyr::select(tmp_data, -sample, -time, -Status)
    n_chans <- length(unique(data$electrode))
  }

  # check for existing reference.
  if (!is.null(data$reference)) {
    tmp_data <- tmp_data + data$reference$ref_data
  }

  # Convert ref_chan channel numbers into channel names
  if (is.numeric(ref_chans)) {
    ref_chans <- names(tmp_data)[ref_chans]
  }

  # Get excluded channel names and/or convert to numbers if necessary
  if (!is.null(exclude)){
    if (is.numeric(exclude)) {
      excluded <- names(tmp_data)[exclude]
    } else {
        excluded <- exclude
        exclude <- which(names(tmp_data) == exclude)
    }
  } else {
    excluded <- NULL
  }

  if (length(ref_chans) == 1 && ref_chans == "average") {
    if (is.null(exclude)) {
      if (robust) {
        ref_data <- matrixStats::rowMedians(as.matrix(tmp_data))
      } else {
        ref_data <- rowMeans(tmp_data)
      }
    } else {
      if (robust) {
        ref_data <- matrixStats::rowMedians(as.matrix(tmp_data[-exclude]))
      } else {
        ref_data <- rowMeans(tmp_data[-exclude])
      }
    }
    tmp_data <- tmp_data - ref_data

  } else {
    if (any(all(ref_chans %in% colnames(tmp_data)) | is.numeric(ref_chans))) {
      if (length(ref_chans) > 1) {
        ref_data <- rowMeans(tmp_data[, ref_chans])
      } else {
        ref_data <- unlist(tmp_data[, ref_chans])
      }
      tmp_data <- tmp_data - ref_data
    } else {
      stop("Electrode(s) not found.")
    }
  }

  if (is.eeg_data(data)) {
    data$signals <- tmp_data
    if (any(ref_chans == "average")) {
      data$reference <- list(ref_chans = ref_chans, ref_data = ref_data,
                             excluded = excluded)
    } else {
      data$reference <- list(ref_chans = ref_chans, ref_data = ref_data,
                             excluded = NULL)
    }
    return(data)
  } else {
    return(as_tibble(tmp_data))
  }
}

#' Baseline correction.
#'
#' Used to remove the mean of a specified time period from the data. Currently
#' only performs subtractive baseline. With a data frame, searches for
#' "electrode" and "epoch" columns, and groups on these when found. An electrode
#' column is always required; an epoch column is not.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to be baseline corrected.
#' @param time_lim Numeric character vector (e.g. time_lim <- c(-.1, 0)). If
#'   none given, defaults to mean of whole epoch if the data is epoched, or the
#'   channel mean if the data is continuous.
#' @import dplyr
#' @import tidyr
#' @export

rm_baseline <- function(data, time_lim = NULL) {

  obj_class <- is.eeg_data(data)

  if (obj_class) {
    orig_data <- data
    data <- cbind(data$signals, data$timings)
    if (orig_data$continuous){
      data <- tidyr::gather(data, electrode, amplitude, -time,
                            -sample, factor_key = TRUE)
    } else {
      data <- tidyr::gather(data, electrode, amplitude, -time,
                            -epoch, -sample, factor_key = TRUE)
    }
  }

  if (!("time" %in% colnames(data))) {
    stop("Time dimension is required.")
  }

  if (length(time_lim) == 1) {
    stop("time_lim should specify the full time range.")
  }

  # if the data is epoched, group by electrode and epoch; otherwise, just by
  # electrode.
  if ("epoch" %in% colnames(data)) {
    data <- dplyr::group_by(data, electrode, epoch, add = TRUE)
  } else{
    data <- dplyr::group_by(data, electrode, add = TRUE)
  }

  if (is.null(time_lim)) {
    # if no time_lim provided, just delete mean of all time points
    data <- dplyr::mutate(data, amplitude = amplitude - mean(amplitude))
  } else {

    data_sel <- dplyr::filter(data,
                              time >= time_lim[1],
                              time <= time_lim[2])
    baseline <- dplyr::summarise(data_sel,
                                 bl = mean(amplitude))
    # This is relatively memory intensive - not so bad now but would prefer
    # another way. Could get extremely painful with time-frequency data.
    data <- dplyr::left_join(data, baseline)
    data <- dplyr::mutate(data, amplitude = amplitude - bl)
    data <- dplyr::select(data, -bl)
  }

  data <- ungroup(data)

  if (obj_class) {
    data <- tidyr::spread(data, electrode, amplitude)
    if (orig_data$continuous) {
      orig_data$signals <- dplyr::select(data, -time, -sample)
      orig_data$timings <- dplyr::select(data, time, sample)
    } else {
      orig_data$signals <- dplyr::select(data, -time, -epoch, -sample)
      orig_data$timings <- dplyr::select(data, time, epoch, sample)
    }
    orig_data
  } else {
    data
  }
}

#' Create epochs from EEG data
#'
#' Creates epochs around specified event triggers. Requires data of class
#' \code{eeg_data}.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Continuous data to be epoched.
#' @param ... Parameters passed to functions
#'
#' @export

epoch_data <- function(data, ...) {
  UseMethod("epoch_data", data)
}

#' Create epochs from EEG data
#'
#' @param data Continuous data to be epoched.
#' @param ... Parameters passed to functions
#' @export
#'

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

epoch_data.eeg_data <- function(data, events, time_lim = c(-1, 1), ...) {

  if (!all(events %in% unique(data$events$event_type))) {
    stop("Events not found - check event codes.")
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

  samps <- seq(round(time_lim[[1]] * srate),
               round(time_lim[[2]] * (srate - 1)),
               by = samp_diff)

  event_table <- data$events

  epoch_zero <-
    sort(unlist(purrr::map(events,
                           ~ event_table[which(event_table$event_type == .), ]$event_onset)))

  epoched_data <- purrr::map(epoch_zero,
                             ~ . + samps)

  epoched_data <- purrr::map_df(epoched_data,
                                ~ tibble::tibble(sample = ., time = samps / srate),
                                .id = "epoch")

  # create new event_table NOW
  event_table <- dplyr::inner_join(event_table, epoched_data,
                                   by = c("event_onset" = "sample"))

  epoched_data$epoch <- as.numeric(epoched_data$epoch)
  epoched_data <- dplyr::left_join(epoched_data,
                                   cbind(data$signals, data$timings),
                                   by = c("sample" = "sample"))


  if (!is.null(data$reference)) {
    ref_data <- dplyr::left_join(epoched_data,
                                 as.data.frame(
                                   cbind(sample = data$timings$sample,
                                         ref_data = data$reference$ref_data)),
                                 by = c("sample" = "sample"))
    data$reference$ref_data <- ref_data[["ref_data"]]
  }

  epoched_data$time.y <- NULL
  names(epoched_data)[[3]] <- "time"
  data$signals <- epoched_data[, -1:-3]
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
#' @author Matt Craddock \email{matt@mattcraddock.com}
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
  }

  message(paste0("Downsampling from ", data$srate, "Hz to ",
                 data$srate / q, "Hz."))

  # make sure any saved reference gets downsampled too.
  if (!is.null(data$reference)) {
    data$signals["ref_data"] <- data$reference$ref_data
  }

  # step through each column and decimate each channel
  data$signals <- purrr::map_df(as.list(data$signals),
                                ~signal::decimate(., q))

  # separate reference from main data again
  if (!is.null(data$reference)) {
    data$reference$ref_data <- data$signals["ref_data"]
    data$signals["ref_data"] <- NULL
  }

  # select every qth timing point, and divide srate by
  data$srate <- data$srate / q
  data$timings <- data$timings[seq(1, length(data$timings[[1]]), by= q), ]

  # The event table also needs to be adjusted. Note that this inevitably jitters
  # event timings by up to q/2 sampling points.
  nearest_samps <- findInterval(data$events$event_onset,
                                data$timings$sample)
  data$events$event_onset <- data$timings$sample[nearest_samps]
  data$events$event_time <- data$timings$time[nearest_samps]
  data
}


#' Combine EEG objects
#'
#' Combine multiple \code{eeg_epochs} or \code{eeg_data} objects into a single
#' object. Note that this does not currently perform any sort of checking for
#' internal consistency or duplication. It simply combines the objects in the order they
#' are passed.
#'
#' @param data An \code{eeg_epochs} object
#' @param ... addtional \code{eeg_epochs} objects
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export
#'

eeg_combine <- function(data, ...) {
  UseMethod("eeg_combine", data)
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
      if (!is.null(data$reference$ref_data)) {
        data$reference$ref_data <- c(data$reference$ref_data,
                                     sapply(args,
                                            function(x) x$reference$ref_data))
      }
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
    if (!is.null(data$reference$ref_data)) {
      data$reference$ref_data <- c(data$reference$ref_data,
                                   sapply(args,
                                          function(x) x$reference$ref_data))
    }
  } else {
    stop("All inputs must be eeg_epochs objects.")
  }

  data
}
