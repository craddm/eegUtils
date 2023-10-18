#'@rdname ar_FASTER
#'@export
eeg_FASTER <- function(data, ...) {
  .Deprecated("ar_FASTER",
              msg = "Please use ar_FASTER; eeg_FASTER will be removed in a forthcoming release.")
  ar_FASTER(data, ...)
}
