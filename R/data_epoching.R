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
#' @return Returns an epoched object of class \code{eeg_epochs}
#'
#' @describeIn epoch_data Epoch \code{eeg_data} objects
#'
#' @export

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

  if (is.null(data$epochs)) {
    data$epochs$recording <- NA
  }

  epoch_trigs <- event_table[event_table$event_type %in% events, c("event_type", "epoch")]
  epochs <- tibble::tibble(epoch = unique(epoched_data$epoch),
                           recording = as.character(data$epochs$recording))
  epochs <- dplyr::left_join(epochs, epoch_trigs, by = "epoch")
  data <- eeg_epochs(data = epoched_data[, -1:-3],
                     srate = data$srate,
                     timings = epoched_data[, 1:3],
                     events = event_table,
                     chan_info = data$chan_info,
                     reference = data$reference,
                     epochs = epochs)
  #data$signals <- epoched_data[, -1:-3] #check which columns this is...
  #data$timings <- epoched_data[, 1:3]
  #data$continuous <- FALSE
  #data$events <- event_table
  # data$epochs <- data.frame(epochs = unique(data$timings$epoch),
  #                           recording = )
  # class(data) <- c("eeg_epochs", "eeg_data")
  data
}

#' @export

epoch_data.eeg_epochs <- function(data, ...) {
  stop("Data is already epoched, cannot currently re-epoch.")
}
