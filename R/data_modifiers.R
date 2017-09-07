#' Referencing
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
    n_chans <- dim(tmp_data)[[2]]
  } else {
    tmp_data <- tidyr::spread(data, electrode, amplitude)
    tmp_data <- dplyr::select(tmp_data, -sample, -time, -Status)
    n_chans <- length(unique(data$electrode))
  }
  if (ref_chans == "all") {
    tmp_data <- tmp_data - rowMeans(tmp_data)
  } else {
    if (ref_chans %in% colnames(tmp_data)) {

    } else {
      stop("Electrode(s) not found.")
    }
  }
  if (is.eeg_data(data)) {
    data$signals <- tmp_data
    data
  } else {
    tmp_data
  }
}

#' Baseline correction.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data Data to be baseline corrected.
#' @param time_lim Numeric character vector (e.g. time_lim <- c(-.1, 0)). If none given, defaults to mean of whole epoch if the data is epoched, or the channel mean if the data is continuous.
#' @import dplyr
#' @import tidyr
#' @export

rm_baseline <- function(data, time_lim = NULL) {

  obj_class <- is.eeg_data(data)

  if (obj_class) {
    orig_data <- data
    data <- cbind(data$signals, data$timings)
    data <- tidyr::gather(data, electrode, amplitude, -time, -epoch, -samp_no)
  }

  if (!("time" %in% colnames(data))) {
    stop("Time dimension is required.")
  }

  if (length(time_lim) == 1) {
    stop("time_lim should specify the full time range.")
  }

  # if the data is epoched, group by electrode and epoch; otherwise, just by electrode.
  if("epoch" %in% colnames(data)) {
    data <- dplyr::group_by(data, electrode, epoch)
  } else{
    data <- dplyr::group_by(data, electrode)
  }

  if (is.null(time_lim)){ # if no time_lim provided, just delete mean of all time points
    data <- dplyr::mutate(data, amplitude = amplitude - mean(amplitude))
  } else {
    #
    data_sel <- dplyr::filter(data, time >= time_lim[1], time <= time_lim[2])
    baseline <- dplyr::summarise(data_sel, bl = mean(amplitude))
    # This is relatively memory intensive - not so bad now but would prefer another way.
    # Could get extremely painful with time-frequency data.
    data <- dplyr::left_join(data, baseline)
    data <- dplyr::mutate(data, amplitude = amplitude-bl)
    data <- dplyr::select(data, -bl)
  }

  data <- ungroup(data)

  if (obj_class) {
    data <- spread(data, electrode, amplitude)
    orig_data$signals <- select(data, -time, -epoch, -samp_no)
    orig_data$timings <- select(data, time, epoch, samp_no)
    orig_data
  } else {
    data
  }
}

#' Create epochs.
#'
#' Requires data of class \code{eeg_data}.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param data Continous data to be epoched.
#' @param events Character vector of events to epoch around.
#' @param time_lim Time in seconds to form epoch around the events. Defaults to one second either side.
#' @importFrom dplyr left_join
#' @importFrom purrr map map_df
#'
#' @return Returns an epoched object of class \code{eeg_data}
#'
#' @export

epoch_data <- function(data, events, time_lim = (c(-1,1))) {

  if (is.eeg_data(data)) {
    samps <- seq(round(time_lim[[1]] * data$srate), round(time_lim[[2]] * (data$srate - 1)))
    event_table <- data$events
    epoch_zero <- sort(unlist(map(events, ~event_table[which(event_table$event_type == .), ]$event_onset)))
    epoched_data <- purrr::map(epoch_zero, ~. + samps)
    epoched_data <- purrr::map_df(epoched_data,~data.frame(sample = ., time = samps/data$srate), .id = "epoch")
    epoched_data <- dplyr::left_join(epoched_data, cbind(data$signals, data$timings), by = c("sample" = "sample"))
    epoched_data$time.y <- NULL
    names(epoched_data)[[3]] <- "time"
    data$signals <- epoched_data[,-1:-3]
    data$timings <- epoched_data[,1:3]
    data$continuous <- FALSE
    return(data)
  } else {
    stop("Currently only working for objects of class eeg_data.")
  }
}
