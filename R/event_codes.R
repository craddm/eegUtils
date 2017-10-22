#' Tag events
#'
#' Give trigger events meaningful labels.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @param data
#'
#'

tag_events <- function(data) {

}

#' List events
#'
#' List trigger types and any labels found in an \code{eeg_data} object.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data
#'
#'
#'

list_events <- function(data) {
  unique(data$events$event_type)
}
