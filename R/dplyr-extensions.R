#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom dplyr filter
#' @export
filter.eeg_epochs <- function(.data,
                              ...) {

  orig_cols <- channel_names(.data)
  args <- rlang::enexprs(...)

  # convert the signals to a data frame that has both the EEG data and epoch
  # labels etc. filter out anything that matches the criteria then return to
  # original format. May want to recode this?

  .data$signals <- tibble::as_tibble(.data)
  .data$signals <- dplyr::filter(.data$signals,
                                 ...)
  .data$signals <- .data$signals[, orig_cols]

  lhs_args <- unlist(lapply(args,
                            function(x) as.character(x[[2]])))

  if (is.null(.data$epochs)) {
    warning("Epochs structure missing; Update your eeg_epochs object using update_eeg_epochs.")
    return(.data)
  }

  in_epochs <- lhs_args %in% names(epochs(.data))

  if (any(in_epochs)) {
    .data$epochs <- dplyr::filter(.data$epochs,
                                  !!!args[in_epochs])
    keep_epochs <- unique(.data$epochs$epoch)
    .data$timings <- dplyr::filter(.data$timings,
                                   epoch %in% keep_epochs)
    events(.data) <- dplyr::filter(events(.data),
                                   epoch %in% keep_epochs)
  }
  .data

  in_timings <- lhs_args %in% names(.data$timings)

  if (any(in_timings)) {
    .data$timings <- dplyr::filter(.data$timings,
                                   !!!args[in_timings])
    events(.data) <- dplyr::filter(events(.data),
                                   !!!args[in_timings])
  }

  in_events <- lhs_args %in% names(.data$events)

  if (any(in_events)) {
    .data$events <- dplyr::filter(.data$events,
                                  !!!args[in_events])
  }

  .data
}

#' @importFrom dplyr filter
#' @export
filter.eeg_data <- function(.data, ...) {

  args <- rlang::enexprs(...)

  arg_list <- parse_args(args, .data)

  orig_cols <- channel_names(.data)

  .data$signals <- as.data.frame(.data)
  .data$signals <- dplyr::filter(.data$signals, ...)
  .data$signals <- .data$signals[, orig_cols]

  if (any(arg_list$in_timings)) {
    .data$timings <- dplyr::filter(.data$timings, ...)
    .data$events$time <- .data$events$event_time
    .data$events <- dplyr::filter(.data$events,
                                  ...)
  }

  # # ensure this also handles the epoch structure correctly
  # epo_args <- grepl(paste(names(.data$epochs), collapse = "|"),
  #                   unlist(args))

  if (any(arg_list$in_epochs)) {
    .data$epochs <- dplyr::filter(.data$epochs,
                                  !!!args[arg_list$in_epochs])
  }
  .data
}

#' @importFrom dplyr filter
#' @export
filter.eeg_evoked <- function(.data,
                              ...) {

  orig_cols <- channel_names(.data)
  args <- rlang::exprs(...)
  arg_list <- parse_args(args,
                         .data)

  .data$signals <- as.data.frame(.data)

  is_group_df <- inherits(.data,
                          "eeg_group")

  .data$signals <- dplyr::filter(.data$signals,
                                 ...)
  .data$signals <- .data$signals[, orig_cols]

  if (nrow(.data$signals) == 0) {
    warnings("All data removed!")
  }

  if (any(arg_list$in_timings)) {
    .data$timings <- dplyr::filter(.data$timings, ...)
  }

  #conditionally filter the epochs structure if any of the arguments refer to
  #its contents. also need to fix this for timings - if we filter based on the
  #epochs structure, it may miss out the timings structure
  if (is.null(.data$epochs)) {
    warning("Epochs structure missing; Update your eeg_epochs object using update_eeg_epochs.")
    return(.data)
  }

  if (any(arg_list$in_epochs)) {
    .data$epochs <- dplyr::filter(.data$epochs,
                                  !!!args[arg_list$in_epochs])
    if (is_group_df) {
      part_epo_list <- paste0(.data$epochs$participant_id,
                              .data$epochs$epoch)
      timing_list <- paste0(.data$timings$participant_id,
                            .data$timings$epoch)
      .data$timings <- .data$timings[timing_list %in% part_epo_list, ]
    }

  }
  .data
}

#' @export
filter.eeg_tfr <- function(.data, ...) {

  args <- rlang::enexprs(...)
  args_done <- logical(length(args))

  which_calls <- vapply(args,
                        is.call,
                        logical(1),
                        USE.NAMES = FALSE)

  if (!all(which_calls)) {
     error_message <-
       paste("Invalid call - did you use a named argument instead of a relational operator? (e.g. epoch = 10)")
     stop(error_message,
          call. = FALSE)
  }

  lhs_args <- unlist(lapply(args[which_calls],
                              function(x) as.character(x[[2]])))
  mat_dims <- names(dimnames(.data$signals)) %in% lhs_args
  hmz <- dimnames(.data$signals)[mat_dims]
  hmz <- lapply(hmz,
                function(x) {
                  if (any(grepl("[A-Z]", x))) {x}
                  else {as.numeric(x)}
                })

  in_epochs <- lhs_args %in% names(epochs(.data))

  if (any(in_epochs)) {
    epochs(.data) <- dplyr::filter(epochs(.data),
                                   !!!args[in_epochs])
    keep_epochs <- epochs(.data)$epoch
    .data$timings <- dplyr::filter(.data$timings,
                                   epoch %in% keep_epochs)
    .data$signals <- abind::asub(.data$signals,
                                 as.character(keep_epochs),
                                 dims = which(.data$dimensions == "epoch"),
                                 drop = FALSE)
    .data$events <- dplyr::filter(.data$events,
                                  epoch %in% keep_epochs)
    args_done[in_epochs] <- TRUE
  }

  in_timings <- lhs_args %in% names(.data$timings)

  if (any(in_timings)) {
    .data$timings <- dplyr::filter(.data$timings,
                                   !!!args[in_timings])
    keep_times <- unique(.data$timings$time)
    time_idx <- which(hmz$time %in% unique(keep_times))
    .data$signals <- abind::asub(.data$signals,
                                 time_idx,
                                 dims = which(.data$dimensions == "time"),
                                 drop = FALSE)
    args_done[in_timings] <- TRUE
  }

  in_frequency <- lhs_args == "frequency"

  if (any(in_frequency)) {
    frequency <- as.numeric(dimnames(.data$signals)$frequency)
    args_done[in_frequency] <- TRUE
    freq_args <- which(in_frequency)
    logi_freq <- lapply(freq_args,
                        function(x) rlang::eval_tidy(args[[x]]))

    logi_freq <- Reduce("&", logi_freq)
    .data$signals <- abind::asub(.data$signals,
                                 logi_freq,
                                 dims = which(names(dimnames(.data$signals)) %in% "frequency"),
                                 drop = FALSE)
    .data$freq_info$freqs <- .data$freq_info$freqs[logi_freq]
    #.data$signals[[]]
  }

  if (any(args_done == FALSE)) {
    warning(paste("Some arguments not used:",
                  unname(vapply(args[!args_done],
                                rlang::quo_text,
                                character(1)))))
  }
  .data
}

#' This function checks which elements of the object need to be
#' accessed/filtered
#' @keywords internal
parse_args <- function(arg_list,
                       data) {

  # which_calls <- vapply(arg_list,
  #                       is.call,
  #                       logical(1),
  #                       USE.NAMES = FALSE)
  #
  # lhs_args <- unlist(lapply(arg_list[which_calls],
  #                           function(x) as.character(x[[2]])))

  lhs_args <- unlist(
    lapply(arg_list,
           all.vars)
    )

  in_epochs <- lhs_args %in% names(epochs(data))
  in_timings <- lhs_args %in% names(data$timings)

  return(list(lhs_args = lhs_args,
              in_epochs = in_epochs,
              in_timings = in_timings))
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
  .data$signals <- dplyr::select(.data$signals,
                                 ...)
  .data$mixing_matrix <- dplyr::select(.data$mixing_matrix,
                                       ...,
                                       electrode)
  keep_comps <- channel_names(.data)
  .data$unmixing_matrix <- dplyr::filter(.data$unmixing_matrix,
                                         Component %in% keep_comps)

  if (!is.null(.data$chan_info)) {
    .data$chan_info <- .data$chan_info[.data$chan_info$electrode %in% .data$mixing_matrix$electrode, ]
  }
  .data
}

#' @importFrom dplyr select
#' @export
select.eeg_stats <- function(.data, ...) {
  .data$statistic <- dplyr::select(.data$statistic,
                                   ...)
  .data$pvals <- dplyr::select(.data$pvals,
                               ...)
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

#' @importFrom dplyr rename
#' @export
dplyr::rename

#' @importFrom dplyr rename
#' @export
rename.eeg_ICA <- function(.data,
                           ...) {
  .data$signals <- dplyr::rename(.data$signals,
                                 ...)
  .data$mixing_matrix <- dplyr::rename(.data$mixing_matrix,
                                       ...)
  .data$unmixing_matrix$Component <- names(.data$signals)
  .data
}

#' @importFrom dplyr rename
#' @export
rename.eeg_epochs <- function(.data,
                           ...) {
  .data$signals <- dplyr::rename(.data$signals,
                                 ...)
  .data$chan_info$electrode <- names(.data$signals)
  .data
}


