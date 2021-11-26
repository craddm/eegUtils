#' Modify the epochs structure
#'
#' Get or set the epochs structure of an `eegUtils` object.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data `eegUtils` object to view
#' @examples
#'   epochs(demo_spatial)
#' @export
epochs <- function(data) {
  if (any(class(data) %in% c("eeg_data",
                             "eeg_epochs",
                             "eeg_ICA",
                             "eeg_tfr",
                             "eeg_evoked"))) {
    data$epochs
  }
}


#' @param value Structure to replace `epochs` structure with.
#' @rdname epochs
#' @export
`epochs<-` <- function(data,
                       value) {

  if (any(class(data) %in% c("eeg_data",
                             "eeg_epochs",
                             "eeg_ICA",
                             "eeg_tfr",
                             "eeg_evoked"))) {
    data$epochs <- value
  }
  data
}

#' Query and set elements of the `epochs` metadata structures
#'
#' These functions are used to query and set the `participant_id` and
#' `recording` fields of an `epochs` structure within any `eegUtils` object.
#' Note that the two `set_*` functions do not operate on `eeg_group` objects.
#'
#' @param data An `eegUtils` object.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @return * `get_participant_id()` returns a character vector containing the
#' `participant_id`
#' * `get_recording()` returns a character vector containing
#' the `recording` from values of `recording` from the `eegUtils` object
#' * `set_participant_id()` returns the full `eegUtils` object with the
#' `participant_id` field modified in the `epochs` structure
#' * `set_recording()` returns the full `eegUtils`object with the `recording` field modified in the
#' `epochs` structure
#' @examples
#' get_participant_id(demo_epochs)
#' get_recording(demo_epochs)
#' set_participant_id(demo_epochs, "002")
#' set_recording(demo_epochs, "test_recording")
#' @export
get_participant_id <- function(data) {
  unique(data$epochs$participant_id)
}

#' @rdname get_participant_id
#' @export
get_recording <- function(data) {
  unique(data$epochs$recording)
}

#' @param participant_id A character vector giving the name to use for the
#'   `participant_id` for a particular object.
#' @rdname get_participant_id
#' @export
set_participant_id <- function(data,
                               participant_id) {
  if (inherits(data, "eeg_group")) {
    stop("Cannot set `participant_id` for `eeg_group` objects, only for objects from individual participants.")
  }
  if (!is.character(participant_id)) {
     stop("`participant_id` must be a character vector")
  }
  data$epochs$participant_id <- participant_id
  data
}

#' @rdname get_participant_id
#' @param recording A character vector giving the name to use for the EEG
#'   recording.
#' @export
set_recording <- function(data,
                          recording) {
  if (inherits(data, "eeg_group")) {
    stop("Cannot set `recording` for `eeg_group` objects, only for objects from individual participants.")
  }
  if (!is.character(recording)) {
    stop("`recording` must be a character vector")
  }
  data$epochs$recording <- recording
  data
}
