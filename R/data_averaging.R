#' Calculate averages (e.g. ERPs) for single datasets
#'
#' This function is used to create an \code{eeg_evoked} object from
#' \code{eeg_epochs}.
#'
#' @param data An \code{eeg_epochs} object.
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

#' @describeIn eeg_average Create evoked data from \code{eeg_epochs}
#' @importFrom tibble tibble
#' @importFrom dplyr left_join group_by_at summarise_at ungroup
#' @param cols Columns from the \code{epochs} structure that the average should
#'   group on. NULL, the default, uses all columns other than the \code{epoch}
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
    if ("participant_id" %in% cols) {
      col_names <- cols
    } else {
      col_names <- c("participant_id", cols)
    }
  } else {
    col_names <- names(data$epochs)
    col_names <- col_names[!(col_names %in% c("epoch"))]
  }

  data$signals <-
    dplyr::group_by_at(data$signals,
                       .vars = vars(time, col_names)) %>%
    dplyr::summarise_at(.vars = vars(elecs),
                        mean) %>%
    dplyr::group_by_at(.vars = col_names) %>%
    dplyr::mutate(epoch = dplyr::group_indices()) %>%
    dplyr::ungroup()

  timings <- data$signals[, c("time", "epoch", col_names)]

  epochs <- dplyr::select(timings,
                          epoch,
                          !!col_names) %>%
    dplyr::distinct()

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

#' @describeIn eeg_average average an eeg_tfr objects over epochs.
#' @export
eeg_average.eeg_evoked <- function(data,
                                   cols = NULL,
                                   ...) {
  message("Data is already averaged.")
  data
}

#' @describeIn eeg_average average an eeg_tfr objects over epochs.
#' @export
eeg_average.eeg_tfr <- function(data,
                                cols = NULL,
                                ...) {
  if (!"epoch" %in% data$dimensions) {
    message("Data is already averaged.")
  } else {
    orig_names <- dimnames(data$signals)
    data <- average_tf(data)
    data$dimensions <- data$dimensions[-which(data$dimensions == "epoch")]
  }
  data
}


#' Internal function for averaging over epochs for eeg_tfr objects.
#' @param data data to average over
#' @keywords internal
average_tf <- function(data) {

  # Need to find a way to make this respect epochs structure...
  orig_dims <- dimnames(data$signals)

  if (data$freq_info$output == "phase") {
    data$signals <- apply(data$signals,
                          c(2, 3, 4),
                          circ_mean)

  } else {
    avg_tf <- array(0, dim = dim(data$signals)[2:4])
    for (iz in 1:dim(data$signals)[3]) {
      for (ij in 1:dim(data$signals)[4]) {
        avg_tf[, iz, ij] <- colMeans(data$signals[ , , iz, ij, drop = FALSE])
      }
    }

    data$signals <- avg_tf
    dimnames(data$signals) <- orig_dims[2:4]
    cols <- c("epoch", "participant_id", "recording")
    data$epochs <- data$epochs[1, colnames(data$epochs) %in% cols]
  }
  data$timings <- dplyr::filter(data$timings, epoch == 1)
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
