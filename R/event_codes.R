#' Tag events
#'
#' Give trigger events meaningful labels. Existing labels will be overwritten.
#' Use hierarchical labelling to tag an event with multiple labels: separate
#' labels with a "/" symbol. (e.g. "cond1" for a trigger that belongs to one
#' condition, "cond1/cond2" for a trigger that could belong to more than one
#' condition).
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data} or \code{eeg_epochs}
#' @param ... Parameters passed to S3 methods
#' @family event handlers
#' @export

tag_events <- function(data, ...) {
  UseMethod("tag_events", data)
}

#' @param trigs Character vector of trigger numbers
#' @param event_label Labels for the events.
#' @importFrom tibble as_tibble
#' @export
#' @describeIn tag_events Tag events in an \code{eeg_data} object

tag_events.eeg_data <- function(data,
                                trigs,
                                event_label,
                                ...) {

  if (length(trigs) != length(event_label)) {
    stop("Trigs and event_label parameters must be the same length.")
  }

  if (!any(trigs %in% unlist(list_events(data)))) {
    stop(paste0("Trigger(s) not found. Check trigger values with list_events()."))
  }

  data$events <- merge(data$events,
                       data.frame(event_type = trigs,
                                  event_label = as.character(event_label),
                                  stringsAsFactors = FALSE))
  data$events <- tibble::as.tibble(data$events)
  data
}

#' @describeIn tag_events Tag events in an epoched dataset
#' @export
tag_events.eeg_epochs <- function(data,
                                  trigs,
                                  event_label,
                                  ...) {

  if (length(trigs) != length(event_label)) {
    stop("Trigs and event_label parameters must be the same length.")
  }

  data$events <- merge(data$events,
                       data.frame(event_type = trigs,
                                  event_label = as.character(event_label),
                                  stringsAsFactors = FALSE))
  data$events <- tibble::as.tibble(data$events)
  data
}

#' List events
#'
#' List trigger types and any labels found in an \code{eeg_data} object.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data An object of class \code{eeg_data}
#' @examples
#' list_events(demo_epochs)
#' @export
#' @family event handlers
#' @seealso \code{\link{tag_events}} and \code{\link{list_epochs}}

list_events <- function(data) {
  if (!is.eeg_data(data)) {
    stop("For eeg_data objects only.")
  }

  if ("event_label" %in% names(data$events)) {
    data$events[!duplicated(data$events$event_type), c("event_type", "event_label")]
  } else {
    data.frame(event_type = unique(data$events$event_type))
  }

}

#' List epochs
#'
#' List trigger types and any labels found in an \code{eeg_epochs} object.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_epochs}
#' @param ... Additional arguments
#' @export
#' @family event handlers
#' @seealso \code{\link{tag_events}} and \code{\link{list_events}}

list_epochs <- function(data, ...) {
  UseMethod("list_epochs", data)
}

#' @describeIn list_epochs List epochs and associated events from \code{eeg_epochs} objects
#' @export
list_epochs.eeg_epochs <- function(data, ...) {
  data$events[, c("epoch", "event_type", "event_label")]
}

#' @describeIn list_epochs List epochs and associated events from \code{eeg_ICA} objects
#' @export
list_epochs.eeg_ICA <- function(data, ...) {
  data$events[, c("epoch", "event_type", "event_label")]
}

#' Modify events structure
#'
#' Get or set the values in the `events` structure of an eegUtils object.
#'
#' @examples
#' events(demo_epochs)
#' events(demo_epochs) <- mutate(events(demo_epochs),
#'  sf = dplyr::case_when(
#'          event_type %% 2 == 0 ~ "HSF",
#'          event_type %% 2 == 1 ~ "LSF",
#'  ))
#' events(demo_epochs)
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param .data \code{eegUtils} object to view
#' @export
events <- function(.data) {
  UseMethod("events", .data)
}

#' @export
events.eeg_data <- function(.data) {
  .data$events
}

#' @export
events.eeg_epochs <- function(.data) {
  .data$events
}


#' @param value Value to replace `events` structure with.
#' @rdname events
#' @export
`events<-` <- function(.data, value) {
  UseMethod("events<-", .data)
}

#' @rdname events
#' @export
`events<-.eeg_epochs` <- function(.data, value) {

  .data$events <- value
  .data
}

#' @rdname events
#' @export
`events<-.eeg_data` <- function(.data, value) {
  .data$events <- value
  .data
}

