#' Calculate averages (e.g. ERPs) for single datasets
#'
#' This function is used to create an \code{eeg_evoked} object from \code{eeg_epochs}.
#'
#' @param data An \code{eeg_epochs} object.
#' @param ... Other arguments passed to the averaging functions
#' @author Matt craddock \email{matt@@mattcraddock.com}
#' @export

eeg_average <- function(data, ...) {
  UseMethod("eeg_average", data)
}

#' @describeIn eeg_average Default method for averaging EEG objects
#' @export

eeg_average.default <- function(data, ...) {
  stop("eeg_epochs or eeg_tfr object required as input.")
}

#' @param cond_label Only pick events that include a given label. Character
#'   vector.
#' @param calc_var Can be used to calculate measures of variability around the
#'   average.
#' @describeIn eeg_average Create evoked data from \code{eeg_epochs}
#' @export
eeg_average.eeg_epochs <- function(data,
                                   cond_label = NULL,
                                   calc_var = NULL, ...) {

  if (is.null(cond_label)) {

    n_times <- length(unique(data$timings$time))
    n_epochs <- length(unique(data$timings$epoch))
    orig_elecs <- names(data$signals)
    data$signals <- as.matrix(data$signals)
    dim(data$signals) <- c(n_times,
                           n_epochs,
                           ncol(data$signals))
    data$signals <- base::colMeans(base::aperm(data$signals, c(2, 1, 3)))
    data$signals <- as.data.frame(data$signals)
    names(data$signals) <- orig_elecs

    if (!is.null(calc_var) && calc_var == "SE") {
      data_sd <- lapply(data$signals,
                        function(x) matrixStats::colSds(as.matrix(x)) / sqrt(nrow(x)))
      data_sd <- as.data.frame(do.call("rbind", data_sd))
      names(data_sd) <- names(data$signals)
      data$var <- data_sd
    }

  } else {

    # Check for presence of labels
    lab_check <- label_check(cond_label,
                             unique(list_epochs(data)$event_label))

    if (!all(lab_check)) {
      stop("Not all labels found. Use list_events to check labels.")
    }

    evoked <- vector("list", length(cond_label))

    for (i in seq_along(cond_label)) {
      tmp_dat <- select_epochs(data,
                               epoch_events = cond_label[[i]])
      tmp_dat$signals <- split(tmp_dat$signals,
                               tmp_dat$timings$time)
      tmp_dat$signals <- lapply(tmp_dat$signals,
                                colMeans)
      evoked[[i]] <- as.data.frame(do.call("rbind",
                                           tmp_dat$signals))
      row.names(evoked[[i]]) <- NULL
    }
    data$signals <- evoked

    if (length(cond_label) > 1) {
      names(data$signals) <- cond_label
    } else {
      data$signals <- data$signals[[1]]
    }
  }
  data$events <- NULL
  data$timings <- unique(data$timings["time"])
  class(data) <- c("eeg_evoked")
  data
}

#' @describeIn eeg_average average an eeg_tfr objects over epochs.
#' @export
eeg_average.eeg_tfr <- function(data,
                                cond_label = NULL,
                                calc_var = NULL, ...) {
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
#' @noRd

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
