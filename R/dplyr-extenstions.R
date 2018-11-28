#' Re-export filter from dplyr
#' @importFrom dplyr filter
#' @name filter
#' @rdname filter
#' @export
NULL

#' @importFrom dplyr select filter
#' @export
filter.eeg_epochs <- function(.data, ...) {
  orig_cols <- names(.data$signals)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data
}

#' @importFrom dplyr select filter
#' @export
filter.eeg_data <- function(.data, ...) {
  orig_cols <- names(.data$signals)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data
}

#' @export
dplyr::select

#' @importFrom dplyr select
#' @export
select.eeg_epochs <- function(.data, ...) {
  .data$signals <- dplyr::select(.data$signals, ...)
  if (!is.null(.data$chan_info)) {
    .data$chan_info <- dplyr::filter(.data$chan_info,
                                     electrode %in% names(.data$signals))
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
#' @noRd

mutate.eeg_data <- function(.data, ...) {
  .data$signals <- dplyr::mutate(.data$signals, ...)
  .data
}

#' @importFrom dplyr mutate
#' @noRd
mutate.eeg_epochs <- function(.data, ...) {
  .data$signals <- dplyr::mutate(.data$signals, ...)
  .data
}

