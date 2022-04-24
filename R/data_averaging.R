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
  stop("eeg_epochs or eeg_tfr object required as input.")
}

#' @describeIn eeg_average Create evoked data from `eeg_epochs`
#' @importFrom tibble tibble
#' @importFrom dplyr left_join group_by_at summarise_at ungroup
#' @param cols Columns from the `epochs` structure that the average should
#'   group on. NULL, the default, uses all columns other than the `epoch`
#'   column.
#' @export
eeg_average.eeg_epochs <- function(data,
                                   cols = NULL,
                                   ...) {

  elecs <- channel_names(data)

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

  data$signals <-
    dplyr::summarise(data$signals,
                     dplyr::across(dplyr::all_of(elecs),
                                   mean))

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


#' @describeIn eeg_average average an `eeg_evoked` object over epochs.
#' @export
eeg_average.eeg_evoked <- function(data,
                                   cols = NULL,
                                   ...) {
  is_group_df <- inherits(data, "eeg_group")

  if (is.null(cols)) {
    message("Data is already averaged.")
    return(data)
  } else {
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
                      dplyr::across(
                        c(time,
                          dplyr::all_of(col_names))
                        )
                      )
    data$signals <-
      dplyr::summarise(data$signals,
                       dplyr::across(dplyr::all_of(elecs),
                                     mean))
    # We want to end up with a unique epoch number for each level of the main
    # variable we are grouping by - so not by participant_id or by time. e.g.
    # each participant should have an epoch 1, and it be the same epoch 1. So if
    # grouping by 2 categories, each combination should have a unique epoch number.
    # Only a concern for group data.
    data$signals <-
      dplyr::group_by(data$signals,
                      dplyr::across(
                        dplyr::all_of(cols)
                        )
                      )
    data$signals <-
      dplyr::mutate(data$signals,
                    epoch = dplyr::cur_group_id())

    data$signals <- dplyr::ungroup(data$signals)
    timings <- data$signals[, c("time", "epoch", "participant_id", col_names)]

    epochs <- dplyr::select(timings,
                            epoch,
                            !!col_names)
    epochs <- unique(epochs)
    if (is_group_df) {
      timings <- data$signals[, c("time", "epoch", "participant_id")]
    } else {
      timings <- data$signals[, c("time", "epoch")]
    }
    timings <- unique(timings)

    # if (!("recording" %in% colnames(epochs))) {
    #   epochs$recording <- recording_id
    # }
    data$timings <- timings
    data$epochs <- epochs
  }

  data
}

#' @describeIn eeg_average average an `eeg_tfr` object over epochs.
#' @export
eeg_average.eeg_tfr <- function(data,
                                cols = NULL,
                                ...) {
  if (!any(c("participant_id", "epoch") %in% data$dimensions)) {
    message("Data is already averaged.")
  } else {
    data <- average_tf(data,
                       cols = cols)
    class(data) <- c("tfr_average",
                     class(data))
  }
  data
}


#' Internal function for averaging over epochs for eeg_tfr objects.
#' @param data data to average over
#' @keywords internal
average_tf <- function(data,
                       cols = NULL) {

  # Need to find a way to make this respect epochs structure...
  orig_dims <- dimnames(data$signals)

  is_group_df <- inherits(data,
                          "eeg_group")
  if ("participant_id" %in% names(orig_dims)) {
    data$signals <- aperm(data$signals,
                          c("participant_id",
                            "epoch",
                            "time",
                            "electrode",
                            "frequency"))
    orig_dims <- dimnames(data$signals)
    #is_group_df <- TRUE
  }

  if (!is.null(cols)) {
    if ("participant_id" %in% cols) {
      col_names <- cols
    } else {
      col_names <- c(
        "participant_id",
        cols
        )
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

  epo_types <- unique(
    epochs(data)[col_names]
    )
  new_epos <- nrow(epo_types)
  n_times <- length(
    dimnames(data$signals)$time
    )

  # There must be a less hacky way of doing this
  epo_nums <-
    lapply(1:new_epos,
           function(x) dplyr::inner_join(
             epochs(data),
             epo_types[x, , drop = FALSE],
             by = col_names)[["epoch"]])

  # convert epoch numbers from epochs() to positions in matrix
  epo_nums <- lapply(epo_nums,
                     function(x) which(orig_dims$epoch %in% x))

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
    epo_types$epoch <- 1:new_epos
    epochs(data) <- epo_types

  } else {
    if (identical(data$freq_info$output,
                  "phase")) {
      data$signals <- apply(
        data$signals,
        c(2, 3, 4),
        circ_mean
        )
    } else if (identical(data$freq_info$output,
                         "power")) {

    # maybe one day try to calm this nested loop gore down
    # but it's pretty quick so hey
    avg_tf <- array(0, dim = c(new_epos,
                               dim(data$signals)[2:4]))

    for (iz in 1:dim(data$signals)[3]) { # electrodes
      for (ij in 1:dim(data$signals)[4]) { # frequencies
        for (ik in 1:new_epos) {
          avg_tf[ik, , iz, ij] <-
            array(
              colMeans(
                data$signals[epo_nums[[ik]], ,
                             iz, ij,
                             drop = FALSE]),
              dim = c(1, n_times,
                      1, 1)
              )
          }
      }
    }

    data$signals <- avg_tf
    new_dims <- orig_dims
    new_dims[["epoch"]] <- as.character(1:new_epos)
    dimnames(data$signals) <- new_dims
    epo_types$epoch <- 1:new_epos
    epochs(data) <- epo_types
    } else {
      stop("Averaging of fourier coefficients not supported.")
    }
  }
  data$timings <-
    tibble::tibble(
      time = rep(
        as.numeric(dimnames(data$signals)[["time"]]),
        new_epos
        ),
      epoch = rep(1:new_epos,
                  each = n_times)
      )
  data
}

#' @export

eeg_average.eeg_group <- function(data,
                                  ...) {
  NextMethod()
  #stop("Not currently supported for `eeg_group` objects.")
}
