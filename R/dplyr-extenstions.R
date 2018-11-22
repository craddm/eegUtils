#' @noRd
filter.eeg_epochs <- function(.data, ...) {
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals$time <- NULL
  .data$signals$epoch <- NULL
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data
}

#' @noRd
filter.eeg_data <- function(.data, ...) {
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals$time <- NULL
  .data$signals$sample <- NULL
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data
}

#' @noRd
select.eeg_epochs <- function(.data, ...) {
  .data$signals <- dplyr::select(.data$signals, ...)
  if (!is.null(.data$chan_info)) {
    .data$chan_info <- .data$chan_info[.data$chan_info$electrode %in% names(.data$signals), ]
  }
  .data
}

#' @noRd
select.eeg_data <- function(.data, ...) {
  .data$signals <- dplyr::select(.data$signals, ...)
  if (!is.null(.data$chan_info)) {
    .data$chan_info <- .data$chan_info[.data$chan_info$electrode %in% names(.data$signals), ]
  }
  .data
}
