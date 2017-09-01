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
    if (ref_chans %in% colnames(tmp_data)){

    } else {
      stop("Electrode(s) not found.")
    }
  }

  tmp_data
}

#' Baseline correction. If no time_lim
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data Data to be baseline corrected.
#' @param time_lim Numeric character vector (e.g. time_lim <- c(-.1, 0)). If none given, defaults to mean of whole epoch if the data is epoched, or the channel mean if the data is continuous.

rm_baseline <- function(data, time_lim = NULL) {

  if (!("time" %in% colnames(data))) {

  }

  if (is.null(time_lim)){
    time_lim <- c(min(data$time), max(data$time))
  }

}

#' Create epochs. Requires data of class \code{eeg_data}.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data
#' @param events Character vector of events to epoch around.
#' @param time_lim Time in seconds to form epoch around the events. Defaults to one second either side.
#' @importFrom dplyr left_join rename
#' @importFrom purrr map map_df

epoch_data <- function(data, events, time_lim = (c(-1,1))) {
  if (is.eeg_data(data)) {
    samps <- seq(round(time_lim[[1]] * data$srate), round(time_lim[[2]] * (data$srate - 1)))
    event_table <- data$events
    epoch_zero <- sort(unlist(map(events, ~event_table[which(event_table$event_type == .), ]$event_onset)))
    epoched_data <- map(epoch_zero, ~. + samps)
    epoched_data <- map_df(epoched_data,~data.frame(samp_no = ., time = samps/data$srate), .id = "epoch_no")
    epoched_data <- dplyr::left_join(epoched_data, cbind(data$signals, data$timings), by = c("samp_no" = "sample"))
    epoched_data$time.y <- NULL
    names(epoched_data)[[3]] <- "time"
    epoched_data
  } else {
    stop("Currently only working for objects of class eeg_data.")
  }
}
