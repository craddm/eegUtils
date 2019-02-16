#'@rdname ar_FASTER
#'@export
eeg_FASTER <- function(.data, ...) {
  .Deprecated("ar_FASTER",
              msg = "Please use ar_FASTER; eeg_FASTER will be removed in a forthcoming release.")
  ar_FASTER(.data, ...)
}

#' Load EEGLAB .set files
#'
#' DEPRECATED - use \code{import_set}
#'
#' EEGLAB .set files are standard Matlab .mat files, but EEGLAB can be set to
#' export either v6.5 or v7.3 format files. Only v6.5 files can be read with
#' this function. v7.3 files (which use HDF5 format) are not currently
#' supported, as they cannot be fully read with existing tools.
#'
#'
#' @param file_name Filename (and path if not in present working directory)
#' @param df_out Defaults to FALSE - outputs an object of class eeg_data. Set to
#'   TRUE for a normal data frame.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @importFrom R.matlab readMat
#' @importFrom dplyr group_by mutate rename
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr is_empty
#' @export

load_set <- function(file_name, df_out = FALSE) {
  .Deprecated("import_set",
              msg = c("load_set is deprecated and has been renamed to import_set. ",
                      "It will be removed in a forthcoming release."))
  file_out <- import_set(file_name, df_out = FALSE)
  file_out
}


#'@rdname eeg_reference
#'@export
reref_eeg <- function(data, ...) {
  .Deprecated("eeg_reference",
              msg = "reref_eeg has been renamed eeg_reference, and will be removed in a future version of eegUtils.")
  eeg_reference(data, ...)
}
