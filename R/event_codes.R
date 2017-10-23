#' Tag events
#'
#' Give trigger events meaningful labels. Existing labels will be overwritten.
#' Use hierarchical labelling to tag an event with multiple labels: separate
#' labels with a "/" symbol. (e.g. "cond1" for a trigger that belongs to one
#' condition, "cond1/cond2" for a trigger that could belong to more than one
#' condition).
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data}
#' @param trigs Character vector of trigger numbers
#' @param event_label Labels for the events.
#' @importFrom dplyr left_join
#' @importFrom tibble tibble
#' @export
#' @seealso \code{\link{list_events}}

tag_events <- function(data, trigs, event_label) {

  if (is.eeg_data(data)) {
    if (length(trigs) != length(event_label)) {
      error("Trigs and event_label parameters must be the same length.")
    }

    if ("event_label" %in% names(data$events)) {
      data$events <- data$events[-3]
    }

    data$events <- dplyr::left_join(data$events, tibble::tibble(event_type = trigs, event_label = event_label), by = "event_type")
    data
  } else {
      error("Object is not of class eeg_data.")
    }
}

#' List events
#'
#' List trigger types and any labels found in an \code{eeg_data} object.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data}
#'
#' @export
#'
#' @seealso \code{\link{tag_events}}

list_events <- function(data) {
  evs <- unique(data$events$event_type)
  if (event_label %in% names(data$events)) {
    data.frame(event_no = evs)
  } else {
    data.frame(event_no = evs, event_label = unique(data$events$event_label))
  }
}
