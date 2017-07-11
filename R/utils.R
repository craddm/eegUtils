#' Get standard electrode locations
#'
#' Joins standard electrode location information from eegUtils internal data.
#'
#' @param data An EEG dataset.
#' @param electrode The column name containing electrode names in data.
#'
#' @import dplyr
#' @export

join_locs <- function(data, ele = "electrode", .drop = FALSE) {
  data[, ele] <- toupper(unlist(data[, ele]))
  electrodeLocs[, ele] <- toupper(unlist(electrodeLocs[, ele]))

  if (.drop) {
    data <- dplyr::left_join(data, electrodeLocs, by = ele)
    data <- filter(data, !(electrode %in% c("NAZ", "LM", "RM", "VEOG", "HEOG")))
  } else {
    data <- dplyr::right_join(data, electrodeLocs, by = ele)
  }

  return(data)
}
