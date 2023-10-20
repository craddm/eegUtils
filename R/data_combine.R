#' Combine `eegUtils` objects
#'
#' Combine multiple `eeg_epochs`, `eeg_data`, or `eeg_evoked` objects into a
#' single object. The function will check the `participant_id` entry in the
#' `epochs` structure of each object to see if the objects come from a single
#' participant or from multiple participants. If the data are from multiple
#' participants, it will create an `eeg_group` object. For individual
#' participants, it will check for duplicate epochs. For most objects, it will
#' concatenate the objects if none are found. However, for `eeg_data` it will
#' instead try to correct the epoch numbers. Check the details below for further
#' advice.
#'
#' # Combining `eeg_data` objects
#'
#' Combining `eeg_data` is mainly intended to be used for combining multiple
#' recordings from a single participant prior to subsequent epoching. Thus,
#' `check_timings` defaults to true, and the function will change the epochs and
#' timing structures of the resulting combined object to be as if it were a
#' single recording. The objects will be combined in the input order, so ensure
#' that the objects are input in chronological order.
#'
#' # Combining `eeg_epochs` objects
#'
#' There are several scenarios where you might wish to combine `eeg_epochs`. For
#' example, a user may have processed continuous data in smaller chunks
#' reflecting short recording blocks before epoching. They then wish to combine
#' these into a single object. In that case, the epoch numbering should reflect
#' chronological ordering and needs to be corrected.
#'
#' If `check_timings == TRUE`, the function will perform several checks before
#' combining objects. First, it will check for duplicate epochs in the `epochs`
#' structure of each object. If each object only has unique epochs, the objects
#' will be combined without correction. Thus, combining across separate
#' recordings or separate participants will not elicit correction. The user
#' should ensure
#'
#' If there are any duplicates (e.g. a participant has more than one epoch
#' numbered one from the same recording), it will then check if there are any
#' missing epochs. If there are, the new trial numbering cannot be automatically
#' determined, so the objects cannot be combined without further manual
#' intervention. If there are no missing epochs, it will then check if there is
#' any decreases in epoch numbers across objects. If there are any, then the
#' epoch numbers and timings for the objects will be adjusted.
#'
#' Alternatively, the user may wish to combine `eeg_epochs` objects from
#' different participants or from entirely different recording sessions of the
#' same participant. In this case, no correction of timings or epoch numbers is
#' desirable. `check_timings == TRUE` should detect this and skip correction,
#' but can be explicitly set to `FALSE`.
#'
#' @param data An `eeg_data`, `eeg_epochs`, or `eeg_evoked` object, or a list of
#'   such objects.
#' @param ... additional `eeg_data` or `eeg_epochs` objects
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @return If all objects have the same `participant_id`, returns an object of
#'   the same class as the original input object. If the objects have different
#'   `participant_id` numbers, an object of both class `eeg_group` and the same
#'   class as the original input object.
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
    "Cannot combine objects of class ",
    class(data)[[1]],
    call. = FALSE
  )
}

#' @describeIn eeg_combine Method for combining lists of `eeg_data` and
#'   `eeg_epochs` objects.
#' @export
eeg_combine.list <- function(data,
                             ...) {

  list_classes <- lapply(data,
                         class)
  # unname the list as the names mess with combining
  data <- unname(data)

  if (any(is.list(unlist(list_classes)))) {
    stop("Cannot handle list of lists. Check class of list items.")
  }

  if (check_classes(list_classes)) {

    return(do.call(eeg_combine,
                   data))
  } else {
    stop("All list elements must be of the same class.")
  }
}

#' @describeIn eeg_combine Method for combining `eeg_data` objects.
#' @param check_timings Check whether sample times / epoch numbers are
#'   continuously ascending; if not, modify so that they are. Useful when, for
#'   example, combining epochs derived from multiple recording blocks. Defaults to TRUE
#' @export

eeg_combine.eeg_data <- function(data,
                                 ...,
                                 check_timings = TRUE) {
  args <- list(...)

  if (length(args) == 0) {
    stop("Nothing to combine.")
  }

  check_participants(data,
                     args)

  if (all(vapply(args,
                 is.eeg_data,
                 logical(1)))) {
    if (any(sapply(args,
                   is.eeg_epochs))) {
      stop("All objects must be unepoched eeg_data objects.")
    } else {
      nsamps <- samples(data)

      if (length(args) > 1) {
        nsamps <- c(nsamps,
                    vapply(args,
                           samples,
                           numeric(1)))
        nsamps <- cumsum(nsamps)
      }

      data$signals <- dplyr::bind_rows(data$signals,
                                       purrr::map_df(args,
                                                     ~ .$signals))
      for (i in 1:length(args)) {
        events(args[[i]]) <-
          dplyr::mutate(
            events(args[[i]]),
            event_onset = event_onset + nsamps[[i]],
            event_time = (event_onset - 1) / data$srate
          )
      }

      data$events  <- dplyr::bind_rows(data$events,
                                       purrr::map_df(args,
                                                     ~ .$events))
      data$timings <- dplyr::bind_rows(data$timings,
                                       purrr::map_df(args,
                                                     ~ .$timings))

      data$chan_info <- unique(
        dplyr::bind_rows(data$chan_info,
                         purrr::map_df(args,
                                       ~.$chan_info))
      )

      message("Taking first dataset's recording name.")
    }
  }

  if (check_timings == TRUE) {
    data <- check_timings(data)
  }

  data
}

#' @describeIn eeg_combine Method for combining `eeg_epochs` objects
#' @export

eeg_combine.eeg_epochs <- function(data,
                                   ...,
                                   check_timings = TRUE) {

  args <- list(...)

  if (length(args) == 0) {
    stop("Nothing to combine.")
  }

  check_participants(data,
                     args)

  participant_time <- get_time_rows(data,
                                    args)

  channel_consistency <- check_channels(data, args)

  if (channel_consistency > 1) {
    warning(paste("Channels are not consistent across all datasets.",
                  "Check if it is correct to combine them."))
  }

  if (check_classes(c(list(data), args))) {
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
    data$chan_info <- unique(
      dplyr::bind_rows(data$chan_info,
                       purrr::map_df(args,
                                     ~.$chan_info))
    )
  } else {
    stop("All inputs must be unaveraged eeg_epochs objects.")
  }

  if (length(unique(epochs(data)$participant_id)) > 1) {
    message("Multiple participant_ids; creating eeg_group...")
    class(data) <- c("eeg_group", class(data))
    p_ids <- rle(as.character(epochs(data)$participant_id))
    p_ids$lengths <- participant_time#p_ids$lengths * length(unique(data$timings$time))
    data$timings$participant_id <- inverse.rle(p_ids)
  } else {
    #fix epoch numbering for combined objects, but only when there are single
    #participants
    if (check_timings) {
      data <- check_timings(data)
    }
  }
  data
}

#' @describeIn eeg_combine Method for combining `eeg_evoked` objects
#' @export
eeg_combine.eeg_evoked <- function(data,
                                   ...) {
  args <- list(...)
  if (length(args) == 0) {
    stop("Nothing to combine.")
  }

  check_participants(data,
                     args)

  participant_time <- get_time_rows(data,
                                    args)

  channel_consistency <- check_channels(data, args)

  if (channel_consistency > 1) {
    warning(paste("Channels are not consistent across all datasets.",
                  "Check if it is correct to combine them."))
  }

  if (check_classes(c(list(data), args))) {

    data$signals <- dplyr::bind_rows(data$signals,
                                     purrr::map_df(args,
                                                   ~.$signals))
    data$timings <- dplyr::bind_rows(data$timings,
                                     purrr::map_df(args,
                                                   ~.$timings))
    data$epochs  <- dplyr::bind_rows(data$epochs,
                                     purrr::map_df(args,
                                                   ~.$epochs))
    data$chan_info <- unique(
      dplyr::bind_rows(data$chan_info,
                       purrr::map_df(args,
                                     ~.$chan_info))
    )
  } else {
    stop("All inputs must be `eeg_evoked` objects.")
  }

  if (length(unique(epochs(data)$participant_id)) > 1) {
    message("Multiple participant IDs, creating eeg_group.")
    class(data) <- c("eeg_group", "eeg_evoked", "eeg_epochs")
    p_ids <- rle(as.character(epochs(data)$participant_id))
    p_ids$lengths <- participant_time #p_ids$lengths * participant_time#length(unique(data$timings$time))
    data$timings$participant_id <- inverse.rle(p_ids)
  }

  data
}

#' @export
eeg_combine.eeg_tfr <- function(data,
                                ...,
                                check_timings = TRUE) {
  args <- list(...)

  if (length(args) == 0) {
    stop("Nothing to combine.")
  }

  if (check_classes(append(list(data), args))) {

    # issue here: only handles `eeg_tfr` when same n_epochs for each object
    # assumes when combining that they are averages from different participants
    # so won't work if combining within participants, for example

    new_dim <- length(dim(data$signals)) + 1
    orig_dims <- dimnames(data$signals)

    if (!check_epochs(append(list(data), args))) {
      stop("Cannot combine `eeg_tfr` when n_epochs is different across objects.
           You may need to average across single-trials first using `eeg_average()`")
    }

    data$signals <- abind::abind(
      data$signals,
      do.call(abind::abind,
              list(purrr::map(
                args, ~ .$signals
                ),
                along = new_dim)),
      along = new_dim,
      use.first.dimnames = TRUE)

    data$epochs  <- dplyr::bind_rows(data$epochs,
                                     purrr::map_df(args,
                                                   ~ .$epochs))

    unique_part_ids <- unique(data$epochs$participant_id)

    if (length(unique_part_ids) > 1) {
      dimnames(data$signals) <- c(orig_dims,
                                  list(participant_id = unique_part_ids))
      data$dimensions <- c(data$dimensions,
                           "participant_id")
      class(data) <- c("eeg_group",
                       "eeg_tfr")
      message("Multiple participant IDs, creating eeg_group.")

    } else {
      if (check_timings == TRUE) {
        data <- check_timings(data)
      }
      dimnames(data$signals) <- c(orig_dims,
                                  list(epoch = data$epochs$epoch))
      data$epoch <- c(data$dimensions,
                      "epoch")
      class(data) <- c("eeg_tfr")
    }

  } else {
    stop("All inputs must be eeg_tfr objects.")
  }
  data
}

#' @export
eeg_combine.tfr_average <- function(data,
                                    ...) {

  # with tfr_averages, the epochs may not be in the same order, or may not
  # be present in all data, so the arrays all need to be reorganised

  args <- list(...)

  if (length(args) == 0) {
    stop("Nothing to combine.")
  }

  if (!check_dims(c(list(data),
                    args))) {
    stop("Signal dimensions are mismatched. All objects must have the same number of channels, frequencies, epochs, and times.")
  }

  all_data <- rearrange_tfr(data,
                            args)

  new_dim <- length(dim(data$signals)) + 1
  orig_dims <- dimnames(data$signals)

  data$signals <- do.call(abind::abind,
                          list(purrr::map(
                            all_data, ~ .$signals
                            ),
                            along = new_dim,
                            use.first.dimnames = TRUE))

  data$epochs  <- purrr::map_df(all_data,
                                ~ .$epochs)

  unique_part_ids <- unique(data$epochs$participant_id)

  if (length(unique_part_ids) > 1) {
    dimnames(data$signals) <- c(orig_dims,
                                list(participant_id = unique_part_ids))
    data$dimensions <- c(data$dimensions,
                         "participant_id")
    class(data) <- c("eeg_group",
                     "eeg_tfr")
    message("Multiple participant IDs, creating eeg_group.")
  }
  data
}

#' @export
eeg_combine.eeg_ICA <- function(data,
                                ...) {
  stop("Combination of eeg_ICA objects is not currently supported.")
}

#' Check consistency of event and timing tables
#'
#' @param data `eeg_data` or `eeg_epochs` object
#' @param verbose Print informative messages
#' @keywords internal

check_timings <- function(.data,
                          verbose = TRUE) {
  UseMethod("check_timings", .data)
}

#' @rdname check_timings
#' @keywords internal
check_timings.eeg_data <- function(.data,
                                   verbose = TRUE) {
  .data$timings$sample <- 1:nrow(.data$timings)
  .data$timings$time <- (.data$timings$sample - 1) / .data$srate
  .data
}

#' @rdname check_timings
#' @keywords internal
check_timings.eeg_epochs <- function(data,
                                     verbose = TRUE) {

  n_rows <- nrow(data$timings)

  # Check for duplicate epochs. If not, there's nothing to correct.
  if (!any(duplicated(data$epochs[, c("epoch", "participant_id", "recording")]))) {
    if (verbose) message("No duplicate epochs found, combining objects.")
    return(data)
  }

  # If there are any positive differences other than 1, then some trials have
  # been removed and the trials are not necessarily in chronological order. It
  # then isn't clear how to combine them.
  if (any(diff(data$timings$epoch) > 1)) {
    stop("Some epochs appear to be missing. ",
         "Cannot automatically combine objects. ",
         "See help for further details.")
  }

  # if the epoch numbers are not ascending, fix them...
  if (any(diff(data$timings$epoch) < 0)) {
    if (verbose) {
      message("Replacing timings with ascending values.",
              "This assumes that the objects are provided in the correct order.")
    }
  }

  # Do I need while here, or just if? double check...
  while (any(diff(data$timings$epoch) < 0)) {

    # check consistency of the data timings table.
    # timings should be consistently increasing.
    # only works correctly with 2 objects

    #check for any places where epoch numbers decrease instead of increase
    switch_locs <- which(diff(data$timings$epoch) < 0)

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

#' Rearrange and combine `tfr_average` objects
#'
#' Epochs in different `tfr_average` objects are not necessarily in the same order.
#'
#' @param data first object to be combined
#' @param add_obj list of additional objects
#' @keywords internal
rearrange_tfr <- function(data,
                          add_obj) {

  # Combine all epochs() metadata from the objects to be combined with the first
  # data
  epo_list <- purrr::map_df(add_obj,
                            epochs)

  # Merge them into a single data frame that can be used for subsequent sorting
  sort_list <- merge(epochs(data),
                     epo_list,
                     all = TRUE)

  epo_cols <- colnames(sort_list)

  sort_list <- split(sort_list,
                     sort_list$participant_id,
                     drop = TRUE)
  new_orders <- lapply(sort_list,
                       `[[`, "epoch")
  full_list <- c(list(data),
                 add_obj)
  # Step through list, sorting each one accordingly
  sorted_list <-
    lapply(seq_along(full_list),
           function(x) {
             epochs(full_list[[x]]) <- sort_list[[x]]
             epochs(full_list[[x]])$epoch <- 1:nrow(epochs(full_list[[x]]))
             full_list[[x]]$signals <- full_list[[x]]$signals[new_orders[[x]], , , , drop = FALSE]
             full_list[[x]]
           })
  sorted_list
}

check_dims <- function(x) {

  sig_dims <- lapply(x,
                     function(y) dim(y$signals))
  length(unique(sig_dims)) == 1
  #add more informative error messages
}


check_participants <- function(data,
                               args) {

  part_ids <- c(get_participant_id(data),
                sapply(args, get_participant_id))

  if (any(is.na(part_ids))) {
    stop(
      "`participant_id` is missing from at least one object. ",
      "Please ensure that the `participant_id` field is set in all `eeg_epochs` objects. ",
      "See ?set_participant_id for assistance."
    )
  } else if (any(part_ids == "")) {
    warning(
      "`participant_id` is set to a default value of \"\" in at least one object. ",
      "If your data is from multiple participants, ",
      "please ensure that data from each participant has a unique ",
      "`participant_id` field before using `eeg_combine`. ",
      "See ?set_participant_id for assistance."
    )
  }
}

get_time_rows <- function(data,
                          args) {
  c(nrow(data$timings),
    vapply(args,
           function(x) nrow(x$timings),
           numeric(1))
    )
}

check_channels <- function(data,
                           args) {
  length(
    unique(
      c(list(channel_names(data)),
        lapply(args, FUN = channel_names)))
  ) == 1
}
