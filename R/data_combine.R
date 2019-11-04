#' Combine EEG objects
#'
#' Combine multiple \code{eeg_epochs} or \code{eeg_data} objects into a single
#' object. Note that this does not currently perform any sort of checking for
#' internal consistency or duplication. It simply combines the objects in the
#' order they are passed.
#'
#' @param data An \code{eeg_data} or \code{eeg_epochs} object, or a list of such
#'   objects.
#' @param ... additional \code{eeg_data} or \code{eeg_epochs} objects
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @importFrom dplyr mutate bind_rows
#' @importFrom purrr map_df
#' @export
#'
eeg_combine <- function(data,
                        ...) {
  UseMethod("eeg_combine", data)
}

#'@export
eeg_combine.default <- function(data,
                                ...) {
  stop(
    "Don't know how to combine objects of class",
    class(data)[[1]],
    call. = FALSE
  )
}

#' @describeIn eeg_combine Method for combining lists of \code{eeg_data} and
#'   \code{eeg_epochs} objects.
#' @export
eeg_combine.list <- function(data,
                             ...) {

  out_dat <- data[[1]]
  out_class <- class(data[[1]])
  out_dat$signals <- purrr::map_df(data,
                                   ~.$signals)
  out_dat$events  <- purrr::map_df(data,
                                   ~.$events)
  out_dat$timings <- purrr::map_df(data,
                                   ~.$timings)
  out_dat$epochs <- purrr::map_df(data,
                                  ~.$epochs)
  if (length(unique(epochs(out_dat)$participant_id)) > 1) {
    class(out_dat) <- c("eeg_group", out_class)
  }
  out_dat
}

#' @describeIn eeg_combine Method for combining \code{eeg_data} objects.
#' @export

eeg_combine.eeg_data <- function(data,
                                 ...){

  args <- list(...)
  if (length(args) == 0) {
    stop("Nothing to combine.")
  }
  if (all(vapply(args, is.eeg_data, logical(1)))) {
    if (any(sapply(args, is.eeg_epochs))) {
      stop("All objects must be unepoched eeg_data objects.")
    } else {

      nsamps <- samples(data)

      if (length(args) > 1) {

        nsamps <- c(nsamps,
                    vapply(args,
                           samples,
                           numeric(1)))
      }

      data$signals <- dplyr::bind_rows(data$signals,
                                       purrr::map_df(args,
                                                     ~.$signals))
      for (i in 1:length(args)) {
        events(args[[i]]) <-
          dplyr::mutate(events(args[[i]]),
                        event_onset = event_onset + nsamps[[i]],
                        event_time = (event_onset - 1) / data$srate)
        }

      data$events  <- dplyr::bind_rows(data$events,
                                       purrr::map_df(args,
                                                     ~.$events))
      data$timings <- dplyr::bind_rows(data$timings,
                                       purrr::map_df(args,
                                                     ~.$timings))
      message("Taking first dataset's recording name.")
    }
  }
  data <- check_timings(data)
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
                                     purrr::map_df(args,
                                                   ~.$signals))
    data$events  <- dplyr::bind_rows(data$events,
                                     purrr::map_df(args,
                                                   ~.$events))
    data$timings <- dplyr::bind_rows(data$timings,
                                     purrr::map_df(args,
                                                   ~.$timings))
    data$epochs  <- dplyr::bind_rows(data$epochs,
                                     purrr::map_df(args,
                                                   ~.$epochs))
  } else {
    stop("All inputs must be eeg_epochs objects.")
  }
  #fix epoch numbering for combined objects
  data <- check_timings(data)
  if (length(unique(epochs(data)$participant_id)) > 1) {
    class(data) <- c("eeg_group", class(data))
  }
  data
}


eeg_combine.eeg_evoked <- function(data,
                                   ...) {
  args <- list(...)
  if (length(args) == 0) {
    stop("Nothing to combine.")
  }

  if (check_classes(args)) {

    data$signals <- dplyr::bind_rows(data$signals,
                                     purrr::map_df(args,
                                                   ~.$signals))
    data$timings <- dplyr::bind_rows(data$timings,
                                     purrr::map_df(args,
                                                   ~.$timings))
    data$epochs  <- dplyr::bind_rows(data$epochs,
                                     purrr::map_df(args,
                                                   ~.$epochs))
  } else {
    stop("All inputs must be eeg_evoked objects.")
  }
  class(data) <- c("eeg_GA", "eeg_evoked", "eeg_epochs")
  data
}

#' Check consistency of event and timing tables
#'
#' @param data \code{eeg_data} or \code{eeg_epochs} object
#' @keywords internal

check_timings <- function(.data) {
  UseMethod("check_timings", .data)
}


#' @rdname check_timings
#' @keywords internal
check_timings.eeg_data <- function(.data) {
  .data$timings$sample <- 1:nrow(.data$timings)
  .data$timings$time <- (.data$timings$sample - 1) / .data$srate
  .data
}

#' @rdname check_timings
#' @keywords internal
check_timings.eeg_epochs <- function(data) {

  n_rows <- nrow(data$timings)
  epochs <- unique(data$timings$epoch)

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
      locs <- (switch_locs + 1):n_rows
      data$timings$epoch[locs] <- data$timings$epoch[locs] + switch_epo
      data$timings$sample[locs] <- data$timings$sample[locs] + switch_sample
    } else {
      for (i in 1:(length(switch_locs) - 1)) {
        switch_epo <- data$timings$epoch[switch_locs[i]]
        switch_sample <- data$timings$sample[switch_locs[i]]
        locs <- (switch_locs[i] + 1):switch_locs[i + 1]
        data$timings$epoch[locs] <- data$timings$epoch[locs] + switch_epo
        data$timings$sample[locs] <- data$timings$sample[locs] + switch_sample
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
                                    by = c("epoch",
                                           "time"))
    data$events$event_onset <- data$events$sample
    data$events$sample <- NULL
  }
  data$epochs$epoch <- unique(data$events$epoch)
  data
}
