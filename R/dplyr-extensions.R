#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom dplyr filter
#' @export
filter.eeg_epochs <- function(.data,
                              ...) {

  orig_cols <- channel_names(.data)
  args <- rlang::exprs(...)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals,
                                 ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings,
                                 ...)
  .data$events <- dplyr::filter(.data$events,
                                ...)

  #conditionally filter the epochs structure if any of the arguments refer to
  #its contents. also need to fix this for timings - if we filter based on the
  #epochs structure, it may miss out the timings structure
  if (is.null(.data$epochs)) {
    warning("Epochs structure missing; Update your eeg_epochs object using update_eeg_epochs.")
    return(.data)
  }

  epo_args <- grepl(paste(names(.data$epochs), collapse = "|"),
                    unlist(args))
  if (any(epo_args)) {
    .data$epochs <- dplyr::filter(.data$epochs,
                                  !!!args[epo_args])
  }
  .data
}

#' @importFrom dplyr filter
#' @export
filter.eeg_data <- function(.data, ...) {
  orig_cols <- channel_names(.data)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings, ...)
  .data$events <- dplyr::filter(.data$events,
                                ...)

  # ensure this also handles the epoch structure correctly
  epo_args <- grepl(paste(names(.data$epochs), collapse = "|"),
                    unlist(args))
  if (any(epo_args)) {
    .data$epochs <- dplyr::filter(.data$epochs,
                                  !!!args[epo_args])
  }
  .data
}

#' @importFrom dplyr filter
#' @export
filter.eeg_evoked <- function(.data,
                              ...) {

  orig_cols <- channel_names(.data)
  args <- rlang::exprs(...)
  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals,
                                 ...)
  .data$signals <- .data$signals[, orig_cols]
  .data$timings <- dplyr::filter(.data$timings,
                                 ...)
  .data$events <- dplyr::filter(.data$events,
                                ...)

  #conditionally filter the epochs structure if any of the arguments refer to
  #its contents. also need to fix this for timings - if we filter based on the
  #epochs structure, it may miss out the timings structure
  if (is.null(.data$epochs)) {
    warning("Epochs structure missing; Update your eeg_epochs object using update_eeg_epochs.")
    return(.data)
  }

  epo_args <- grepl(paste(names(.data$epochs), collapse = "|"),
                    unlist(args))
  if (any(epo_args)) {
    .data$epochs <- dplyr::filter(.data$epochs,
                                  !!!args[epo_args])
  }
  .data
}


#' @importFrom dplyr select
#' @export
dplyr::select

#' @importFrom dplyr select
#' @export
select.eeg_epochs <- function(.data,
                              ...) {
  .data$signals <- dplyr::select(.data$signals,
                                 ...)
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

#' @importFrom dplyr select filter
#' @export
select.eeg_ICA <- function(.data, ...) {
  .data$signals <- dplyr::select(.data$signals, ...)
  .data$mixing_matrix <- dplyr::select(.data$mixing_matrix, ..., electrode)
  unmix <- data.table::transpose(.data$unmixing_matrix[, 1:(ncol(.data$unmixing_matrix)-1)])
  names(unmix) <- .data$unmixing_matrix$Component
  unmix <- dplyr::select(unmix, ...)
  remaining_comps <- names(unmix)
  unmix <- data.table::transpose(unmix)
  unmix$Component <- remaining_comps
  names(unmix) <- names(.data$unmixing_matrix)
  .data$unmixing_matrix <- unmix

  if (!is.null(.data$chan_info)) {
    .data$chan_info <- .data$chan_info[.data$chan_info$electrode %in% .data$mixing_matrix$electrode, ]
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
