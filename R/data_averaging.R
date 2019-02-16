#' Calculate averages (e.g. ERPs) for single datasets
#'
#' This function is used to create an \code{eeg_evoked} object from \code{eeg_epochs}.
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
#' @export
eeg_average.eeg_epochs <- function(data,
                                   ...) {

  elecs <- channel_names(data)

  if (is.null(data$epochs)) {
    data$epochs <- tibble::tibble(epoch = unique(data$timings$epoch),
                                  recording = NA)
  }

  data$signals <- dplyr::left_join(cbind(data$signals,
                                         data$timings),
                                   data$epochs, by = "epoch")
  col_names <- names(data$epochs)
  col_names <- col_names[!(col_names %in% c("epoch", "recording"))]

  data$signals <-
    dplyr::group_by_at(data$signals,
                       .vars = vars(time, col_names)) %>%
    dplyr::summarise_at(.vars = vars(elecs),
                        mean) %>%
    dplyr::ungroup()

  data <-
    eeg_evoked(data = data$signals[, elecs],
               chan_info = data$chan_info,
               srate = data$srate,
               timings = data$signals[, c("time", col_names)],
               epochs = NULL)
  data
}

#' @describeIn eeg_average average an eeg_tfr objects over epochs.
#' @export
eeg_average.eeg_tfr <- function(data,
                                ...) {
  if (!"epoch" %in% data$dimensions) {
    #message("Data is already averaged.")
  } else {
    data <- average_tf(data)
    data$dimensions <- data$dimensions[-which(data$dimensions == "epoch")]
  }
  data
}

#' Grand average
#'
#' @param data A list of objects to be averaged over; currently only supports
#'   lists of \code{eeg_evoked} objects
#' @param keep_indivs Keep averages for individual participants. Logical.
#'   Defaults to TRUE.
#' @noRd

eeg_grandaverage <- function(data,
                             keep_indivs = TRUE) {

  if (!class(data) == "list") {
    stop("A list of objects is required.")
  }

  if (!is.eeg_evoked(data[[1]])) {
    stop("Only eeg_evoked objects are supported at this time.")
  }

  # Check for consistency of input
  if (check_classes(data)) {
    if (is.null(dim(data[[1]]$signals))) {
      if (check_conds(data)) {
        conds <- names(data[[1]]$signals)
        ga_sigs <- lapply(seq_along(conds),
                          function(x) create_grandavg(data,
                                                      keep_indivs,
                                                      x))
        names(ga_sigs) <- conds
        grand_avg <- eeg_GA(ga_sigs,
                            srate = data[[1]]$srate,
                            timings = data[[1]]$timings,
                            chan_info = data[[1]]$chan_info,
                            indivs = keep_indivs)
        } else {
          stop("Some conditions are not present in every object.")
          }
      } else {
        ga_sigs <- create_grandavg(data, keep_indivs)
        }
    } else {
      stop("Some objects are of different classes.")
      }
  grand_avg <- eeg_GA(ga_sigs,
                      srate = data[[1]]$srate,
                      timings = data[[1]]$timings,
                      chan_info = data[[1]]$chan_info,
                      indivs = keep_indivs)
  grand_avg
}


#' Create a grand average file from multiple evoked objects
#' @importFrom dplyr bind_rows
#' @keywords internal
create_grandavg <- function(data,
                            keep_indivs,
                            x = NULL) {

  if (is.null(x)) {
    indiv_data <- lapply(data,
                         function(i) i$signals)
  } else {
    indiv_data <- lapply(data,
                         function(i) i$signals[[x]])
  }
  if (keep_indivs) {
    indiv_data <- dplyr::bind_rows(indiv_data,
                                   .id = "sub_no")
    } else {
      indiv_data <- Reduce("+",
                           indiv_data) / length(indiv_data)
    }
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
