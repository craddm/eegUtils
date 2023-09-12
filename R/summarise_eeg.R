#' Calculate simple summary statistics for `eeg_*` objects
#'
#' Calculate the timepoint-by-timepoint mean, standard deviation, standard
#' error, or variance `eeg_epochs` objects.
#'
#' @param data An `eegUtils` object.
#' @param ... Various arguments passed to specific functions
#' @return A tibble
#' @export
eeg_summarise <- function(data,
                          ...) {
  UseMethod("eeg_summarise", data)
}

#' @export
eeg_summarise.default <- function(data,
                                  ...) {
  stop("Not implemented for objects of class ",
       paste(class(data), collapse = "/"))
}

#' @param statistic The statistic to calculate at each timepoint. Defaults to
#'   "sem"
#' @param conditions Conditions to group the data by.
#' @param time_lim Timepoint(s) to summarise. Can be a range, for which a
#'   summary statistic will be provided for each timepoint, or a list of
#'   individual times. If none is supplied, the function will calculate a
#'   summary for every timepoint.
#' @describeIn eeg_summarise Calculate summary statistics for `eeg_epochs`
#'   objects
#' @export
eeg_summarise.eeg_epochs <- function(data,
                                     statistic = c("sem",
                                                   "mean",
                                                   "sd",
                                                   "var"),
                                     conditions = NULL,
                                     time_lim = NULL,
                                     ...) {

  if (!is.null(conditions)) {
    groups <- c(conditions, "time")
  } else {
    groups <- "time"
  }

  chans <- channel_names(data)
  data <- as.data.frame(data)

  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim)
  }

  data <- dplyr::group_by(data,
                          dplyr::across({{ groups }}))

  summary_function <-
    switch(statistic,
           mean = function(x) mean(x),
           sem = function(x) sd(x) / length(x),
           sd = function(x) sd(x),
           var = function(x) var(x)
           )

  data <- dplyr::summarise(
    data,
    dplyr::across(chans,
                  summary_function,
                  .names = paste(statistic,
                                 sep = "_",
                                 "{.col}")
                  )
    )
  data
}

#' Calculate the standard error of the mean
#'
#' @param data An `eegUtils` object
#' @param time_lim A vector of time-points to calculate SEM over.
#' @param conditions Split the calculation across multiple conditions
#' @keywords internal

calculate_sem <- function(data,
                          time_lim = NULL,
                          conditions = NULL) {

  conditions <- rlang::enquo(conditions)
  chans <- channel_names(data)

  if (is.list(time_lim)) {
    all_sems <- lapply(time_lim,
                       function(x) calculate_sem(data,
                                                 time_lim = x,
                                                 conditions = {{conditions}}))
    return(
      do.call("rbind",
              all_sems)
    )

  } else {
    data <- select_times(data,
                         time_lim)
  }
  data <- as.data.frame(data)

  if (!is.null(conditions)) {
    data <- group_by(data,
                     {{conditions}},
                     .data$epoch)
  }
  data <- dplyr::summarise(
    data,
    dplyr::across(chans,
                  mean))
  data <- dplyr::summarise(
    data,
    count = length(unique(.data$epoch)),
    dplyr::across(chans,
                  list(sem = ~sd(.) / sqrt(count))))
  data$time_lim <- list(time_lim)
  data
}

calculate_sd <- function(data,
                         time_lim,
                         conditions = NULL) {
  conditions <- rlang::enquo(conditions)
  chans <- channel_names(data)
  data <- select_times(data,
                       time_lim)
  data <- as.data.frame(data)
  if (!is.null(conditions)) {
    data <- group_by(data,
                     {{conditions}})
  }
  data <- dplyr::summarise(
    data,
    count = length(unique(.data$epoch)),
    dplyr::across(chans,
                  list(sd = ~sd(.))))
  data$time_lim <- list(time_lim)
  data
}
