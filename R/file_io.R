#' Function for reading raw data.
#'
#' Currently only BDF or EDF files are supported, and this only provides a wrapper around the \code{edfReader} package.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @param file_name File to import.
#' @import edfReader
#' @import tools
#' @importFrom purrr map_df


import_raw <- function(file_name, file_path = NULL) {
  file_type <- file_ext(file_name)
  if (file_type == "bdf" | file_type == "edf") {
    data <- readEdfSignals(readEdfHeader(file_name))
    sigs <- map_df(data, "signals")
    srate <- data[[1]]$sRate
    data <- eeg_data(sigs, srate)
  } else {
    warning("Unsupported filetype")
    return()
  }
  return(data)
}
