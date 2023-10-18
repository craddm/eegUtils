#' Calculate averages (e.g. event-related potentials) for single datasets
#'
#' This function is used to create an `eeg_evoked` object from `eeg_epochs`. By
#' default, it will try to keep different conditions in the data separate using
#' the `epochs` metadata from the object, thus yielding one average per
#' condition. Alternatively, the user can specify which averages they want using
#' the `cols` argument.
#'
#' @param data An `eeg_epochs` of `eeg_tfr` object.
#' @param ... Other arguments passed to the averaging functions
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @returns An object of class `eeg_evoked` if applied to `eeg_epochs`;
#'   `eeg_tfr` if applied to `eeg_tfr`.
#' @export

eeg_average <- function(data,
                        ...) {
  UseMethod("eeg_average", data)
}

#' @describeIn eeg_average Default method for averaging EEG objects
#' @export

eeg_average.default <- function(data,
                                ...) {
  stop(paste("Cannot average an object of class", class(data)))
}

#' @describeIn eeg_average Create evoked data from `eeg_epochs`objects
#' @param cols Columns from the `epochs` structure that the average should
#'   group on. NULL, the default, uses all columns other than the `epoch`
#'   column.
#' @param verbose Print informative messages during averaging. Defaults to TRUE
#' @examples
#' eeg_average(demo_spatial)
#' eeg_average(demo_spatial, cols = "everything")
#' @export
eeg_average.eeg_epochs <- function(data,
                                   cols = NULL,
                                   verbose = TRUE,
                                   ...) {

  elecs <- channel_names(data)

  # this is here for compatibility - early versions of the structure have no
  # $epochs entry
  if (is.null(data$epochs)) {
    n_epochs <- length(unique(data$timings$epoch))
    data$epochs <-
      tibble::new_tibble(
        list(
          epoch = unique(data$timings$epoch),
          participant_id = character(n_epochs),
          recording = character(n_epochs)
        ),
        nrow = n_epochs,
        class = "epoch_info"
      )
  }

  data$signals <- dplyr::left_join(cbind(data$signals,
                                         data$timings),
                                   data$epochs, by = "epoch")

  recording_id <- unique(data$signals$recording)

  if (length(recording_id) > 1) {
    recording_id <- recording_id[[1]]
  }

  if (!is.null(cols)) {
    if (identical(cols, "everything")) {
      col_names <- "participant_id"
    } else {
      if ("participant_id" %in% cols) {
        col_names <- cols
      } else {
        col_names <- c("participant_id", cols)
      }
    }
  } else {
    col_names <- names(data$epochs)
    col_names <-
      col_names[!(col_names %in% c("epoch", "recording", "event_type"))]
  }

  # break down into individual calls using updated syntax
  data$signals <-
    dplyr::group_by(data$signals,
                    dplyr::across(c(time,
                                    dplyr::all_of(col_names))))

  if (verbose) {
    message("Creating epochs based on combinations of variables: ",
            paste(col_names, ""))
  }

  # Add epoch weights
  data$signals <-
    dplyr::summarise(data$signals,
                     dplyr::across(dplyr::all_of(elecs),
                                   mean),
                     weight = dplyr::n())

  data$signals <-
    dplyr::group_by(data$signals,
                    dplyr::across(dplyr::all_of(col_names)))

  data$signals <-
    dplyr::mutate(data$signals,
                  epoch = dplyr::cur_group_id())

  data$signals <- dplyr::ungroup(data$signals)

  timings <- data$signals[, c("time", "epoch", col_names)]

  epochs <- dplyr::select(timings,
                          epoch,
                          !!col_names)
  epochs$weight <- data$signals$weight
  epochs <- unique(epochs)

  timings <- data$signals[, c("time", "epoch")]
  timings <- unique(timings)

  if (!("recording" %in% colnames(epochs))) {
    epochs$recording <- recording_id
  }

  class(epochs) <- c("epoch_info",
                     "tbl_df",
                     "tbl",
                     "data.frame")

  data <-
    eeg_evoked(
      data = data$signals[, elecs],
      chan_info = data$chan_info,
      srate = data$srate,
      timings = timings,
      reference = data$reference,
      epochs = epochs
    )
  data
}


#' @param weighted Produce a weighted average over epochs, which accounts for
#'   upstream differences in the number of epochs that contribute to each
#'   average.
#' @param verbose Print informative messages during averaging. Defaults to TRUE
#' @describeIn eeg_average average an `eeg_evoked` object over epochs.
#' @export
eeg_average.eeg_evoked <- function(data,
                                   cols = NULL,
                                   weighted = TRUE,
                                   verbose = TRUE,
                                   ...) {

  is_group_df <- inherits(data,
                          "eeg_group")

  if (is.null(cols) && verbose) {
    message("Data is already averaged - you must specify columns to group your averages by.")
    return(data)
  }

  if (identical(cols, "everything")) {
    col_names <- "participant_id"
  } else {
    if ("participant_id" %in% cols) {
      col_names <- cols
    } else {
      col_names <- c("participant_id", cols)
    }
  }
  elecs <- channel_names(data)
  data$signals <- as.data.frame(data)
  data$signals <-
    dplyr::group_by(data$signals,
      dplyr::across(c(time,
                      dplyr::all_of(col_names))))

  if (verbose) {
    message("Creating epochs based on combinations of variables: ",
            paste(col_names, ""))
  }

  if (weighted) {
    if ("weight" %in% names(data$signals)) {
      if (verbose) message("Calculating weighted means.")
      full_weights <- dplyr::summarise(data$signals,
                                       weight = sum(weight))
      data$signals <-
        dplyr::summarise(data$signals,
                         dplyr::across(dplyr::all_of(elecs),
                                       ~weighted.mean(.x, w = weight)))
      data$signals$weight <- full_weights$weight
      weighted <- TRUE
    } else {
      if (verbose) message("No weights found, calculating unweighted means.")
      data$signals <-
        dplyr::summarise(data$signals,
                         dplyr::across(dplyr::all_of(elecs),
                                       mean))
      weighted <- FALSE
    }
  } else {
    data$signals <-
      dplyr::summarise(data$signals,
                       dplyr::across(dplyr::all_of(elecs),
                                     mean))
    weighted <- FALSE
  }

  # We want to end up with a unique epoch number for each level of the main
  # variable we are grouping by - so not by participant_id or by time. e.g. each
  # participant should have an epoch 1, and it be the same epoch 1. So if
  # grouping by 2 categories, each combination should have a unique epoch
  # number. Only a concern for group data.

  data$signals <-
    dplyr::group_by(data$signals,
                    dplyr::across(
                      dplyr::all_of(col_names)
                    )
    )

  data$signals <-
    dplyr::mutate(data$signals,
                  epoch = dplyr::cur_group_id())

  data$signals <- dplyr::ungroup(data$signals)

  epochs <- dplyr::select(data$signals,
                          epoch,
                          !!col_names)

  if (weighted) epochs$weight <- full_weights$weight

  epochs <- unique(epochs)

  if (is_group_df) {
    timings <- data$signals[, c("time", "epoch", "participant_id")]
  } else {
    timings <- data$signals[, c("time", "epoch")]
  }
  timings <- unique(timings)

  data$timings <- timings
  data$epochs <- epochs
  data$signals <- data$signals[, !colnames(data$signals) %in% c("time", "epoch", "weight", col_names)]
  data
}

#' @describeIn eeg_average average an `eeg_tfr` object over epochs.
#' @export
eeg_average.eeg_tfr <- function(data,
                                cols = NULL,
                                weighted = TRUE,
                                verbose = TRUE,
                                ...) {
  if (!any(c("participant_id", "epoch") %in% data$dimensions)) {
    message("Data is already averaged.")
  } else {
    data <- average_tf(data,
                       cols = cols,
                       weighted = weighted,
                       verbose = verbose)
    class(data) <- c("tfr_average",
                     class(data))
  }
  data
}

#' Internal function for averaging over epochs for `eeg_tfr` objects.
#' @param data data to average over
#' @param weighted Calculate weighted means if TRUE (if possible!)
#' @keywords internal
average_tf <- function(data,
                       cols = NULL,
                       weighted,
                       verbose) {

  # Need to find a way to make this respect epochs structure...
  orig_dims <- dimnames(data$signals)
  is_group_df <- inherits(data, "eeg_group")

  if ("participant_id" %in% names(orig_dims)) {
    data$signals <- aperm(data$signals,
                          c("participant_id",
                            "epoch",
                            "time",
                            "electrode",
                            "frequency"))
    orig_dims <- dimnames(data$signals)
  }

  if (!is.null(cols)) {
    if (identical(cols, "everything")) {
      col_names <- "participant_id"
    } else {
      if ("participant_id" %in% cols) {
        col_names <- cols
      } else {
        col_names <- c("participant_id", cols)
      }
    }
  } else {
    col_names <- names(data$epochs)
    if (is_group_df) {
      col_names <- col_names[!(col_names %in% c("participant_id",
                                                "epoch",
                                                "recording",
                                                "event_type"))]
    } else {
      col_names <- col_names[!(col_names %in% c("epoch",
                                                "recording",
                                                "event_type"))]
    }
  }

  epo_types <- unique(epochs(data)[col_names])
  new_epochs <- nrow(epo_types)
  n_times <- length(dimnames(data$signals)$time)

  # There must be a less hacky way of doing this
  epo_nums <-
    lapply(1:new_epochs,
      function(x) {
        sort(merge(epochs(data),
                   epo_types[x, , drop = FALSE],
                   by = col_names)[, "epoch"])
      }
    )
  # convert epoch numbers from epochs() to positions in matrix
  epo_nums <- lapply(epo_nums,
                     function(x) which(orig_dims$epoch %in% x))

  if ("weight" %in% colnames(epochs(data))) {
    epo_weights <- epochs(data)[["weight"]]
  } else {
    epo_weights <- vapply(epo_nums,
                          length,
                          integer(1))
  }

  if (is_group_df) {
    orig_dims <- dimnames(data$signals)
    orig_dims[["participant_id"]] <- "grand_average"
    data$signals <- colMeans(data$signals)
    dim(data$signals) <- c(1,
                           dim(data$signals))
    dimnames(data$signals) <- orig_dims
    data$dimensions <- c("participant_id",
                         "epoch",
                         "time",
                         "electrode",
                         "frequency")
    epo_types$epoch <- 1:new_epochs
    epochs(data) <- epo_types
  } else {
    if (identical(data$freq_info$output, "phase")) {
      data$signals <- apply(
        data$signals,
        c(2, 3, 4),
        circ_mean
      )
    } else if (identical(data$freq_info$output, "power")) {
      avg_tf <- array(0,
                      dim = c(new_epochs,
                              dim(data$signals)[2:4]))
      final_weights <- integer(new_epochs)

      # figuring out how to do weighted means.

      if (inherits(data, "tfr_average") && weighted) {
        if (verbose) message("Calculating weighted means...")
        if (new_epochs == 1) {
          relative_weights <- epo_weights / sum(epo_weights)
          weighted_means <- sweep(data$signals, c(1,3,4),
                                  relative_weights, "*")
          avg_tf[1, , , ] <- colSums(weighted_means)
          final_weights <- sum(epo_weights)
        } else {
          for (epoch_no in seq_len(new_epochs)) {
            cond_epochs <- epo_nums[[epoch_no]]
            temp <- data$signals[cond_epochs, , , ]
            relative_weights <- epo_weights[cond_epochs] / sum(epo_weights[cond_epochs])
            weighted_means <- sweep(temp, c(1,3,4),
                                    relative_weights, "*")
            avg_tf[epoch_no, , , ] <- colSums(weighted_means)
            final_weights[epoch_no] <- sum(epo_weights[cond_epochs])
          }
        }
      } else {
        for (elec in seq_len(dim(data$signals)[3])) { # electrodes
          for (freq in seq_len(dim(data$signals)[4])) { # frequencies
            for (ik in seq_len(new_epochs)) {
              avg_tf[ik, , elec, freq] <-
                array(
                  colMeans(
                    data$signals[epo_nums[[ik]], ,
                                 elec, freq,
                                 drop = FALSE]),
                  dim = c(1, n_times, 1, 1)
                )
            }
          }
        }
        final_weights <- epo_weights
      }
      data$signals <- avg_tf
      new_dims <- orig_dims
      new_dims[["epoch"]] <- as.character(1:new_epochs)
      dimnames(data$signals) <- new_dims
      epo_types$epoch <- 1:new_epochs
      epochs(data) <- epo_types
      epochs(data)$weight <- final_weights
    } else {
      stop("Averaging of fourier coefficients not supported.")
    }
  }
  data$timings <-
    tibble::tibble(
      time = rep(
        as.numeric(dimnames(data$signals)[["time"]]),
        new_epochs
      ),
      epoch = rep(1:new_epochs,
                  each = n_times)
    )
  data
}

#' @export
eeg_average.eeg_group <- function(data,
                                  ...) {
  NextMethod(...)
  #stop("Not currently supported for `eeg_group` objects.")
}
