#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom dplyr filter
#' @export
filter.eeg_epochs <- function(.data, ...) {
  orig_cols <- names(.data$signals)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data
}

#' @importFrom dplyr filter
#' @export
filter.eeg_data <- function(.data, ...) {
  orig_cols <- names(.data$signals)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data
}

#' @importFrom dplyr select
#' @export
dplyr::select

#' @importFrom dplyr select
#' @export
select.eeg_epochs <- function(.data, ...) {
  .data$signals <- dplyr::select(.data$signals, ...)
  new_cols <- names(.data$signals) # dplyr::filter can't find .data$signals to get names directly
  if (!is.null(.data$chan_info)) {
    .data$chan_info <- dplyr::filter(.data$chan_info,
                                     electrode %in% new_cols)
  }
  .data
}

#' @importFrom dplyr select
#' @export
select.eeg_data <- function(.data, ...) {
  .data$signals <- dplyr::select(.data$signals, ...)
  if (!is.null(.data$chan_info)) {
    .data$chan_info <- .data$chan_info[.data$chan_info$electrode %in% names(.data$signals), ]
  }
  .data
}

#' @importFrom dplyr mutate
#' @export
dplyr::mutate

#' @importFrom dplyr mutate
#' @export

mutate.eeg_data <- function(.data, ...) {
  .data$signals <- dplyr::mutate(.data$signals, ...)
  .data
}

#' @importFrom dplyr mutate
#' @export
mutate.eeg_epochs <- function(.data, ...) {
  .data$signals <- dplyr::mutate(.data$signals, ...)
  .data
}


