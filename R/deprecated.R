#'@rdname ar_FASTER
#'@export
eeg_FASTER <- function(data, ...) {
  .Deprecated("ar_FASTER",
              msg = "Please use ar_FASTER; eeg_FASTER will be removed in a forthcoming release.")
  ar_FASTER(data, ...)
}


#'@rdname eeg_reference
#'@export
reref_eeg <- function(data, ...) {
  .Deprecated("eeg_reference",
              msg = "reref_eeg has been renamed eeg_reference, and will be removed in a future version of eegUtils.")
  eeg_reference(data, ...)
}
