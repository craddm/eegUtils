#' Calculate averages (e.g. ERPs) for single datasets
#'
#' This function is used to create an `eeg_evoked` object from
#' `eeg_epochs`.
#'
#' @param data An `eeg_epochs` object.
#' @param ... Other arguments passed to the averaging functions
#' @author Matt craddock \email{matt@@mattcraddock.com}
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
      tibble::new_tibble(list(epoch = unique(data$timings$epoch),
                              participant_id = character(n_epochs),
                              recording = character(n_epochs)),
                         nrow = n_epochs,
                         class = "epoch_info")
  }

  data$signals <- dplyr::left_join(cbind(data$signals,
                                         data$timings),
                                   data$epochs, by = "epoch")
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
    col_names <- col_names[!(col_names %in% c("epoch", "recording", "event_type"))]
  }

  data$signals <-
    dplyr::group_by_at(data$signals,
                       .vars = vars(time, col_names)) %>%
    dplyr::summarise_at(.vars = vars(elecs),
                        mean) %>%
    dplyr::group_by_at(.vars = col_names) %>%
    dplyr::mutate(epoch = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  timings <- data$signals[, c("time", "epoch", col_names)]

  epochs <- dplyr::select(timings,
                          epoch,
                          !!col_names)
  epochs <- unique(epochs)
  timings <- data$signals[, c("time", "epoch")]
  timings <- unique(timings)

  class(epochs) <- c("epoch_info",
                     "tbl_df",
                     "tbl",
                     "data.frame")

  data <-
    eeg_evoked(data = data$signals[, elecs],
               chan_info = data$chan_info,
               srate = data$srate,
               timings = timings,
               epochs = epochs)
  data
}

#' @describeIn eeg_average average an `eeg_epochs` object over epochs.
#' @export
eeg_average.eeg_evoked <- function(data,
                                   cols = NULL,
                                   ...) {
  message("Data is already averaged.")
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

  if ("participant_id" %in% names(orig_dims)) {
    data$signals <- aperm(data$signals, c("participant_id",
                                          "time",
                                          "electrode",
                                          "frequency"))
    orig_dims <- dimnames(data$signals)
  }


  if (!is.null(cols)) {
    if ("participant_id" %in% cols) {
      col_names <- cols
    } else {
      col_names <- c("participant_id", cols)
    }
  } else {
    col_names <- names(data$epochs)
    col_names <- col_names[!(col_names %in% c("epoch",
                                              "recording",
                                              "event_type"))]
  }

  epo_types <- unique(epochs(data)[col_names])
  new_epos <- nrow(epo_types)
  n_times <- dim(data$signals)[2]

  # There must be a less hacky way of doing this
  epo_nums <-
    lapply(1:new_epos,
           function(x) dplyr::inner_join(epochs(data),
                                         epo_types[x, ],
                                         by = col_names)[["epoch"]])

  # convert epoch numbers from epochs() to positions in matrix
  epo_nums <- lapply(epo_nums,
                     function(x) which(orig_dims$epoch %in% x))

  if (identical(data$freq_info$output, "phase")) {
    data$signals <- apply(data$signals,
                          c(2, 3, 4),
                          circ_mean)

  } else if (identical(data$freq_info$output, "power")) {

    # maybe one day try to calm this nested loop gore down
    # but it's pretty quick so hey
    avg_tf <- array(0, dim = c(new_epos,
                               dim(data$signals)[2:4]))

    for (iz in 1:dim(data$signals)[3]) { # electrodes
      for (ij in 1:dim(data$signals)[4]) { # frequencies
        for (ik in 1:new_epos) {
          avg_tf[ik, , iz, ij] <-
            array(colMeans(data$signals[epo_nums[[ik]], ,
                                        iz, ij,
                                        drop = FALSE]),
                  dim = c(1, n_times,
                          1, 1))
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
  data$timings <-
    tibble::tibble(
      time = rep(as.numeric(dimnames(data$signals)[["time"]]), new_epos),
      epoch = rep(1:new_epos, each = n_times)
      )#dplyr::filter(data$timings, epoch == 1)
  data
}

#' Check that all classes in a list match
#'
#' @param data list of objects to check
#' @keywords internal

check_classes <- function(data) {

  stopifnot(is.list(data))

  dat_classes <- lapply(data,
                        class)
  check_class <- sapply(dat_classes,
                        identical,
                        dat_classes[[1]])
  all(check_class)
}

#' Check that all conditions in a list match
#' @noRd

check_conds <- function(data_list) {

  get_names <- lapply(data_list,
                      function(x) names(x$signals))
  check_names <- sapply(get_names,
                        identical,
                        get_names[[1]])
  all(check_names)
}
