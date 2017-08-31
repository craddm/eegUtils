#' Referencing.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data Data to re-reference. Long format is expected.
#' @param ref_chans Channels to reference data to. Defaults to "all" i.e. average of all electrodes in data.
#' @importFrom tidyr spread
#' @importFrom dplyr select
#'

reref_eeg <- function(data, ref_chans = "all") {
  if (is.eeg_data(data)) {
    tmp_data <- data$signals
    tmp_data <- select(tmp_data, -sample, -time, -Status)
    n_chans <- dim(tmp_data)[[2]]
  } else {
    tmp_data <- spread(data, electrode, amplitude)
    tmp_data <- select(tmp_data, -sample, -time, -Status)
    n_chans <- length(unique(data$electrode))
  }
  if (ref_chans == "all") {
    #non_elecs <- ncol(tmp_data) - n_chans
    #tmp_data[,(non_elecs+1):(dim(tmp_data)[[2]])] <- tmp_data[,-1:-non_elecs ] - rowMeans(tmp_data[,-1:-non_elecs])
    tmp_data <- tmp_data - rowMeans(tmp_data)
  } else {
    if (ref_chans %in% names(tmp_data)){

    } else {
      stop("Electrode(s) not found.")
    }
  }

  tmp_data
}

#' Baseline correction.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data

#' Create epochs.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data
#' @import tidyr
#' @import dplyr
