#' Referencing
#'
#' Used to reference the EEG data to specified electrode or electrodes. Defaults
#' to average reference. Note that this stores the reference channel in the
#' \code{eeg_data} structure; if one already exists, it is added back to the
#' data before further referencing takes place.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Data to re-reference. Primarily meant for use with data of class
#'   \code{eeg_data}. Long format is expected if a data-frame is submitted.
#' @param ref_chans Channels to reference data to. Defaults to "average" i.e.
#'   average of all electrodes in data. Character vector of channel names or numbers.
#' @param exclude Electrodes to exclude from average reference calculation. Only
#' @param robust Use median instead of mean; only applicable for average reference.
#' @importFrom tidyr spread
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom matrixStats rowMedians
#' @return object of class \code{eeg_data}, re-referenced as requested.
#' @export
#'

reref_eeg <- function(data, ref_chans = "average", exclude = NULL, robust = FALSE) {

  if (is.eeg_data(data)) {
    tmp_data <- data$signals
    n_chans <- dim(tmp_data)[[2]]
  } else {
    tmp_data <- tidyr::spread(data, electrode, amplitude)
    tmp_data <- dplyr::select(tmp_data, -sample, -time, -Status)
    n_chans <- length(unique(data$electrode))
  }

  # check for existing reference.
  if (!is.null(data$reference)) {
    tmp_data <- tmp_data + data$reference$ref_data
  }

  # Convert ref_chan channel numbers into channel names
  if (is.numeric(ref_chans)) {
    ref_chans <- names(tmp_data)[ref_chans]
  }

  # Get excluded channel names and/or convert to numbers if necessary
  if (!is.null(exclude)){
    if (is.numeric(exclude)) {
      excluded <- names(tmp_data)[exclude]
    } else {
        excluded <- exclude
        exclude <- which(names(tmp_data) == exclude)
    }
  } else {
    excluded <- NULL
  }

  if (length(ref_chans) == 1 && ref_chans == "average") {
    if (is.null(exclude)) {
      if (robust) {
        ref_data <- matrixStats::rowMedians(as.matrix(tmp_data))
      } else {
        ref_data <- rowMeans(tmp_data)
      }
    } else {
      if (robust) {
        ref_data <- matrixStats::rowMedians(as.matrix(tmp_data[-exclude]))
      } else {
        ref_data <- rowMeans(tmp_data[-exclude])
      }
    }
    tmp_data <- tmp_data - ref_data
  } else {
    if (any(all(ref_chans %in% colnames(tmp_data)) | is.numeric(ref_chans))) {
      if (length(ref_chans) > 1) {
        ref_data <- rowMeans(tmp_data[ , ref_chans])
      } else {
        ref_data <- unlist(tmp_data[ , ref_chans])
      }
      tmp_data <- tmp_data - ref_data
    } else {
      stop("Electrode(s) not found.")
    }
  }

  if (is.eeg_data(data)) {
    data$signals <- tmp_data
    if (any(ref_chans == "average")) {
      data$reference <- list(ref_chans = ref_chans, ref_data = ref_data, excluded = excluded)
    } else {
      data$reference <- list(ref_chans = ref_chans, ref_data = ref_data, excluded = NULL)
    }
    return(data)
  } else {
    return(as_tibble(tmp_data))
  }
}

#' Baseline correction.
#'
#' Used to remove the mean of a specified time period from the data. Currently only performs subtractive baseline.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
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
    if (orig_data$continuous){
      data <- tidyr::gather(data, electrode, amplitude, -time, -sample)
    } else {
      data <- tidyr::gather(data, electrode, amplitude, -time, -epoch, -sample)
    }
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
    data <- tidyr::spread(data, electrode, amplitude)
    if (orig_data$continuous) {
      orig_data$signals <- dplyr::select(data, -time, -sample)
      orig_data$timings <- dplyr::select(data, time, sample)
    } else {
      orig_data$signals <- dplyr::select(data, -time, -epoch, -sample)
      orig_data$timings <- dplyr::select(data, time, epoch, sample)
    }
    orig_data
  } else {
    data
  }
}

#' Create epochs.
#'
#' Creates epochs around specified event triggers. Requires data of class \code{eeg_data}.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param data Continuous data to be epoched.
#' @param events Character vector of events to epoch around.
#' @param time_lim Time in seconds to form epoch around the events. Defaults to one second either side.
#' @importFrom dplyr left_join
#' @importFrom purrr map map_df
#'
#' @return Returns an epoched object of class \code{eeg_data}
#'
#' @export


epoch_data <- function(data, events, time_lim = (c(-1, 1))) {
  if (is.eeg_data(data)) {
    samps <- seq(round(time_lim[[1]] * data$srate), round(time_lim[[2]] * (data$srate - 1)))
    event_table <- data$events
    epoch_zero <- sort(unlist(purrr::map(events,
                                         ~event_table[which(event_table$event_type == .),]$event_onset)))
    epoched_data <- purrr::map(epoch_zero,
                               ~ . + samps)
    epoched_data <-
      purrr::map_df(epoched_data,
                    ~ tibble::tibble(sample = ., time = samps / data$srate),
                    .id = "epoch")
    epoched_data$epoch <- as.numeric(epoched_data$epoch)
    epoched_data <-
      dplyr::left_join(epoched_data,
                       cbind(data$signals, data$timings),
                       by = c("sample" = "sample"))
    if (!is.null(data$reference)) {
      ref_data <-
        dplyr::left_join(epoched_data, as.data.frame(
          cbind(
            sample = data$timings$sample,
            ref_data = data$reference$ref_data
          )
        ), by = c("sample" = "sample"))
      data$reference$ref_data <- ref_data[["ref_data"]]
    }
    epoched_data$time.y <- NULL
    names(epoched_data)[[3]] <- "time"
    data$signals <- epoched_data[, -1:-3]
    data$timings <- epoched_data[, 1:3]
    data$continuous <- FALSE
    return(data)
  } else {
    stop("Currently only working for objects of class eeg_data.")
  }
}
