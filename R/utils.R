#' Get standard electrode locations
#'
#' Joins standard electrode locations to EEG data from eegUtils internal data.
#'
#' @param data An EEG dataset.
#' @param ... Parameters passed to S3 methods.
#' @export

electrode_locations <- function(data, ...) {
  UseMethod("electrode_locations")
}

#' @param electrode The column name containing electrode names in data.
#'   (Defaults to "electrode").
#' @param drop Should electrodes in \code{data} for which default locations are
#'   not available be dropped? (Defaults to FALSE).
#' @param plot Plot obtained electrode locations.
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble is.tibble
#' @describeIn electrode_locations Adds standard locations to a data frame in long format
#' @return A tibble (or data.frame), or ggplot2 object if \code{plot = TRUE}.
#' @export

electrode_locations.data.frame <- function(data,
                                electrode = "electrode",
                                drop = FALSE,
                                plot = FALSE,
                                montage = NULL, ...) {

  if (!is.null(montage)) {
    if (montage == "biosemi64alpha") {
      electrodeLocs[1:64, electrode] <- c(paste0("A", 1:32),
                                          paste0("B", 1:32))
    }
  }

  data[, electrode] <- toupper(data[[electrode]])
  electrodeLocs[, electrode] <- toupper(electrodeLocs[[electrode]])

  if (tibble::is.tibble(data)) {
    elecs <-
      dplyr::pull(unique(data[, electrode])) %in%
      dplyr::pull(electrodeLocs[, electrode])

    if (!all(elecs)) {
      message("Electrodes not found: ",
              paste(dplyr::pull(unique(data[, electrode]))[!elecs], sep = ","))
      } else if (!any(elecs)) {
        stop("No matching electrodes found.")
      }
  } else {
    elecs <-
      unique(data[, electrode]) %in% dplyr::pull(electrodeLocs[, electrode])
    if (!all(elecs)) {
      message("Electrodes not found: ",
              paste(unique(data[, electrode])[!elecs], sep = ","))
    } else if (!any(elecs)) {
      stop("No matching electrodes found.")
    }

  }

  if (drop) {
    data <- dplyr::inner_join(data, electrodeLocs, by = electrode)
  } else {
    data <- dplyr::left_join(data, electrodeLocs, by = electrode)
  }

  if (plot) {
    plotdata <- dplyr::distinct(data, x, y, electrode)
    p <- ggplot2::ggplot(plotdata, aes(x, y)) +
      geom_label(aes(label = electrode))
    return(p)
  } else {
    return(data)
  }
}

#' @param overwrite Overwrite existing channel info. Defaults to FALSE.
#' @import dplyr
#' @import ggplot2
#' @importFrom tibble is.tibble
#' @describeIn electrode_locations Adds standard locations to the chan_info field of an eeg_data object.
#' @export

electrode_locations.eeg_data <- function(data,
                                           drop = FALSE,
                                           plot = FALSE,
                                           montage = NULL,
                                         overwrite = FALSE, ...) {

  if (!is.null(data$chan_info) & !overwrite & !plot) {
    stop("Channel info already present, set overwrite to TRUE to replace.")
  }

  if (!is.null(montage)) {
    if (montage == "biosemi64alpha") {
      electrodeLocs[1:64, "electrode"] <- c(paste0("A", 1:32),
                                          paste0("B", 1:32))
    }
  }

  elec_names <- toupper(names(data$signals))
  electrodeLocs$electrode <- toupper(electrodeLocs$electrode)

  matched_els <- electrodeLocs$electrode %in% elec_names
  missing_els <- elec_names %in% electrodeLocs$electrode

  if (!any(matched_els)) {
    stop("No matching electrodes found.")
  } else if (!all(matched_els)) {
    message(cat("Electrodes not found:", names(data$signals)[!missing_els]))
  }

  data$chan_info <- electrodeLocs[matched_els, ]

  if (drop) {
    data$signals[matched_els]
  }

  if (plot) {
    p <- ggplot2::ggplot(data$chan_info, aes(x, y)) +
      geom_label(aes(label = electrode))
    return(p)
  } else {
    return(data)
  }
}

#' Function to create an S3 object of class \code{eeg_data}
#'
#' @noRd
new_eeg_data <- function(data,
                         srate,
                         events = NULL,
                         chan_info = NULL,
                         timings = NULL,
                         continuous = NULL,
                         reference = NULL) {
  stopifnot(is.double(srate))
  stopifnot(is.data.frame(events))
  structure(signals = data,
            srate = srate,
            events = events)

}

#' Function to create an S3 object of class "eeg_data".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data Raw data - signals from electrodes/channels.
#' @param srate Sampling rate in Hz.
#' @param chan_info String of character names for electrodes.
#' @param timings Timing information - samples and sample /samplirng rate.
#' @param continuous Whether the data is continuous or epoched.
#' @param reference Reference channel information, including names of reference
#'   channels, excluded channels etc.
#' @param events Event table
#' @export

eeg_data <- function(data,
                     srate,
                     events = NULL,
                     chan_info = NULL,
                     timings = NULL,
                     continuous = NULL,
                     reference = NULL) {
  if (srate < 1) {
    stop("Sampling rate must be above 0")
  }
  value <- list(signals = data,
                srate = srate,
                events = events,
                chan_info = chan_info,
                timings = timings,
                continuous = continuous,
                reference = reference
                )
  class(value) <- "eeg_data"
  value
}

#' Function to create an S3 object of class "eeg_evoked"
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data evoked data
#' @param chan_info Electrode locations etc
#' @param timings vector of timepoints
#' @param ... Other parameters

eeg_evoked <- function(data, chan_info, timings, ...) {
  value <- list(signals = data,
                chan_info = chan_info,
                timings = timings)
  class(value) <- c("eeg_evoked", "eeg_data")
}

#' Function to create an S3 object of class "eeg_stats".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param statistic Calculated statistic (e.g. t-statistic)
#' @param pvals calculated p-values for that statistic
#' @param chan_info String of character names for electrodes.
#' @param timings Unique timepoints remaining in the data.
#' @export

eeg_stats <- function(statistic, chan_info, pvals, timings) {

  value <- list(statistic = statistic,
                pvals = pvals,
                chan_info = chan_info,
                timings = timings)
  class(value) <- "eeg_stats"
  value
}

#' Check if object is of class "eeg_data".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object to check.
#' @noRd

is.eeg_data <- function(x) inherits(x, "eeg_data")

#' Check if object is of class "eeg_epochs".
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object to check.
#' @noRd

is.eeg_epochs <- function(x) inherits(x, "eeg_epochs")

#' Check if object is of class \code{eeg_evoked}
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object to check.
#' @noRd
is.eeg_evoked <- function(x) inherits(x, "eeg_evoked")

#' Check if object is of class \code{eeg_stats}
#'
#' @param x Object to check.
#' @noRd
is.eeg_stats <- function(x) inherits(x, "eeg_stats")

#' Check if object is of class \code{eeg_ICA}
#'
#' @param x Object to check.
#' @noRd

is.eeg_ICA <- function(x) inherits(x, "eeg_ICA")

#' Convert eeg_data to data.frame
#'
#' Convert an object of class \code{eeg_data} into a standard data.frame / tibble
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object of class \code{eeg_data}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param events Include events in output.
#' @param ... arguments for other as.data.frame commands

#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_data <- function(x, row.names = NULL,
                                   optional = FALSE, long = FALSE,
                                   events = FALSE, ...) {
  df <- data.frame(x$signals, x$timings)

  if (long) {
    if (x$continuous) {
      df <- tidyr::gather(df,
                          electrode,
                          amplitude,
                          -time,
                          -sample,
                          factor_key = T)
    } else {
      df <- tidyr::gather(df,
                          electrode,
                          amplitude,
                          -time,
                          -sample,
                          -epoch,
                          factor_key = T)
    }
  }

  if (events) {
    df <- dplyr::left_join(df, x$events, by = c("sample" = "event_onset"))
  }
  return(df)
}

#' Convert \code{eeg_epochs} object to data.frame
#'
#' Convert an \code{eeg_epochs} object to a data.frame for use with whatever
#' packages you desire.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_epochs}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param events Include events in output. Defaults to FALSE.
#' @param cond_labels Add column tagging epochs with events that have matching
#'   labels.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_epochs <- function(x, row.names = NULL,
                                     optional = FALSE,
                                     long = FALSE, events = FALSE,
                                     cond_labels = NULL, ...) {

  if (!is.null(cond_labels)) {
    lab_check <- label_check(cond_labels,
                             unique(list_epochs(x))$event_label)
    if (!all(lab_check)) {
      stop("Not all labels found. Use list_events to check labels.")
    }

    df <- lapply(seq_along(cond_labels), function(ix) {
      out_df <- as.data.frame(select_epochs(x, cond_labels[[ix]]))
      out_df$cond_label <- cond_labels[[ix]]
      out_df
    })

    df <- do.call("rbind", df)

  } else {
    df <- data.frame(x$signals,
                     time = x$timings$time,
                     epoch = x$timings$epoch)
    }

    if (long) {
    df <- tidyr::gather(df,
                        electrode,
                        amplitude,
                        -time,
                        -epoch,
                        factor_key = T)
  }

  if (events) {
    df <- dplyr::left_join(df, x$events, by = c("sample" = "event_onset"))
  }
  return(df)
}

#' Convert \code{eeg_evoked} object to data frame
#
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_evoked}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_evoked <- function(x, row.names = NULL,
                                     optional = FALSE, long = FALSE, ...) {
  df <- data.frame(x$signals, time = x$timings$time)
  if (long) {
    df <- tidyr::gather(df,
                        electrode,
                        amplitude,
                        -time,
                        factor_key = T)
  }
  return(df)
}

#' Convert \code{eeg_ICA} object to data frame
#
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object of class \code{eeg_ICA}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_ICA <- function(x, row.names = NULL,
                                  optional = FALSE, long = FALSE, ...) {
  names(x$comp_activations) <- 1:ncol(x$comp_activations)
  df <- data.frame(x$comp_activations, x$timings)
  if (long) {
    df <- tidyr::gather(df,
                        electrode,
                        amplitude,
                        -time,
                        -epoch,
                        -sample,
                        factor_key = T)
  }
  return(df)
}

#' Convert polar to spherical coordinates
#'
#' Hard-coded to a radius of 85 mm (as in BESA).
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param theta Azimuth from polar co-ordinates (theta in supplied electrode
#'   locations)
#' @param phi Elevation from polar co-ordinates (radius in supplied electrode
#'   locations)
#' @noRd

pol_to_sph <- function(theta, phi) {

  theta <- -theta
  phi <- (0.5 - phi) * 180
  sph_theta <- theta / 180 * pi
  sph_phi <- phi / 180 * pi
  radius <- 85

  z <- radius * sin(sph_phi)
  z_cos <- radius * cos(sph_phi)
  x <- z_cos * cos(sph_theta)
  y <- z_cos * sin(sph_theta)
  data.frame(x, y, z)
}

#' Check consistency of events
#'
#' Internal function for checking 1) whether the labels submitted are a mixture
#' of hierarchical and non-hierarchical types 2) whether the labels submitted
#' are present in the data
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param cond_labs labels submitted by the user
#' @param data_labs labels from the actual data
#' @noRd

label_check <- function(cond_labs, data_labs) {

  if (all(grepl("/", cond_labs))) {
    lab_check <- cond_labs %in% data_labs
    } else if (any(grepl("/", cond_labs))) {
      stop("Do not mix hierarchical and non-hierarchical event labels.")
      } else {
        # Check if there is a hierarchical separator "/". If so,
        # split the labels
        if (any(grepl("/", data_labs))) {
          split_labels <- strsplit(data_labs, "/")

          lab_check <- lapply(cond_labs,
                              function(x) vapply(split_labels,
                                                 function(i) x %in% i,
                                                 logical(1)))
          #condense to a single TRUE or FALSE for each label
          lab_check <- vapply(lab_check, any, logical(1))
        } else {
          lab_check <- cond_labs %in% data_labs
        }
      }
}

