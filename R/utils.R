#' Convert eeg_data to data.frame
#'
#' Convert an object of class \code{eeg_data} into a standard data.frame.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_data}
#' @param row.names Kept for compatibility with S3 generic, ignored.
#' @param optional Kept for compatibility with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param events Include events in output.
#' @param coords Include electrode coordinates in output. Only possible when
#'   long = TRUE.
#' @param ... arguments for other as.data.frame commands
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @export

as.data.frame.eeg_data <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   long = FALSE,
                                   events = FALSE,
                                   coords = FALSE,
                                   ...) {
  df <- data.frame(x$signals,
                   x$timings)

  if (long) {
    df <- tidyr::gather(df,
                        electrode,
                        amplitude,
                        -time,
                        -sample,
                        factor_key = T)
    }

  if (events) {
    df <- dplyr::left_join(df,
                           x$events,
                           by = c("sample" = "event_onset"))
  }

  df
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
#' @param events Include events in output. Defaults to FALSE. Currently ignored.
#' @param cond_labels Add column tagging epochs with events that have matching
#'   labels. Deprecated. Metainfo from the epochs structure is now added
#'   automatically.
#' @param coords Include electrode coordinates in output. Ignored if long == FALSE.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @export

as.data.frame.eeg_epochs <- function(x, row.names = NULL,
                                     optional = FALSE,
                                     long = FALSE,
                                     events = FALSE,
                                     cond_labels,
                                     coords = TRUE,
                                     ...) {

  if (!missing(cond_labels)) {
    stop("The cond_labels argument is deprecated.")
  }

  df <- data.frame(x$signals,
                   time = x$timings$time,
                   epoch = x$timings$epoch,
                   stringsAsFactors = FALSE)

  # combine the new data frame with any condition labels from the events table
  if (long) {

    # When converting to long, use factor_key to keep channels in same order,
    # then convert back to character.
    df <- tidyr::gather(df,
                        electrode,
                        amplitude,
                        channel_names(x),
                        factor_key = TRUE)

    df$electrode <- as.character(df$electrode)

    if (coords && !is.null(channels(x))) {
      df <- dplyr::left_join(df,
                             channels(x)[, c("electrode", "x", "y")],
                             by = "electrode")
    }
  }

  if (!is.null(x$epochs)) {
    df <- dplyr::left_join(df,
                           x$epochs,
                           by = "epoch")
  }

  df
}

#' Convert \code{eeg_evoked} object to data frame
#
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_evoked}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param coords include electrode coordinates in output. Ignored if long = FALSE.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_evoked <- function(x,
                                     row.names = NULL,
                                     optional = FALSE,
                                     long = FALSE,
                                     coords = TRUE,
                                     ...) {

  df <- cbind(x$signals,
              x$timings)

  if (long) {
    df <- tidyr::gather(df,
                        "electrode",
                        "amplitude",
                        channel_names(x),
                        factor_key = TRUE)
    df$electrode <- as.character(df$electrode)

    if (coords && !is.null(channels(x))) {
      df <- left_join(df,
                      channels(x)[, c("electrode", "x", "y")],
                      by = "electrode")
    }
  }

  df
}

#' Convert \code{eeg_ICA} object to data frame
#
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_ICA}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param cond_labels add condition labels to data frame. Deprecated.
#' @param mixing If TRUE, outputs the mixing matrix. If FALSE, outputs source activations.
#' @param coords Adds electrode coordinates if TRUE; only if long data and the mixing matrix are requested.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_ICA <- function(x,
                                  row.names = NULL,
                                  optional = FALSE,
                                  long = FALSE,
                                  cond_labels,
                                  mixing = FALSE,
                                  coords = TRUE,
                                  ...) {

  if (!missing(cond_labels)) {
    stop("The cond_labels argument is deprecated.")
  }

  if (mixing) {
    df <- x$mixing_matrix

    if (coords) {

      if (is.null(channels(x))) {
        stop("No chan_info, use electrode_locations() first.")
      }

      df <- dplyr::left_join(df,
                             channels(x)[, c("electrode", "x", "y")],
                             by = "electrode")
    }

  } else {
    df <- data.frame(x$signals,
                     x$timings)

    if (!is.null(x$epochs)) {
      df <- dplyr::left_join(df,
                             x$epochs,
                             by = "epoch")
    }
  }

  if (long) {
    df <- tidyr::gather(df,
                        component,
                        amplitude,
                        channel_names(x),
                        factor_key = T)
    df$component <- as.character(df$component)
  }


  df
}

#' Convert \code{eeg_tfr} objects to data frames
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_tfr}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr spread
#' @export
as.data.frame.eeg_tfr <- function(x,
                                  row.names = NULL,
                                  optional = FALSE,
                                  long = FALSE,
                                  ...) {

  # out_df <- as.data.frame.table(x$signals,
  #                               stringsAsFactors = FALSE)
  if (length(dim(x$signals)) == 3) {
    out_df <- aperm(x$signals, c("time", "frequency", "electrode"))
    dim_size <- dim(out_df)
    dim(out_df) <- c(dim_size[1] * dim_size[2],
                     dim_size[3])
    out_df <- as.data.frame(out_df)
    names(out_df) <- dimnames(x$signals)$electrode
    out_df$time <- rep(as.numeric(dimnames(x$signals)$time),
                       dim_size[2])
    out_df$frequency <- rep(as.numeric(dimnames(x$signals)$frequency),
                            each = dim_size[1])
    out_df$epoch <- 1
  } else {
    out_df <- aperm(x$signals, c("time", "frequency", "epoch", "electrode"))
    dim_size <- dim(out_df)
    dim(out_df) <- c(dim_size[1] * dim_size[2] * dim_size[3],
                     dim_size[4])
    out_df <- as.data.frame(out_df)
    names(out_df) <- dimnames(x$signals)$electrode
    out_df$time <- rep(as.numeric(dimnames(x$signals)$time),
                       dim_size[3] * dim_size[2])#out_df$time)
    out_df$frequency <- rep(as.numeric(dimnames(x$signals)$frequency),
                            each = dim_size[1])
    out_df$epoch <- rep(unique(x$epochs$epoch), each =
                        dim_size[1] * dim_size[2])

    if (!is.null(x$epochs) && "epoch" %in% x$dimensions) {
      #out_df$epoch <- as.numeric(out_df$epoch)
      out_df <- dplyr::left_join(out_df,
                                 x$epochs,
                                 by = "epoch")
    }
  }

  if (long) {
    return(tidyr::gather(out_df, electrode, power, channel_names(x)))
  }

  out_df
}

#' Convert \code{eeg_stats} objects to data frames
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_stats}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param coords Include electrode coordinates in output (ignored if long = FALSE)
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr spread
#' @export

as.data.frame.eeg_stats <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   long = FALSE,
                                   coords = FALSE,
                                   ...) {

  df <- data.frame(x$statistic,
                   time = x$timings)

  if (long) {
    df <- tidyr::gather(df,
                        electrode,
                        statistic,
                        -time,
                        factor_key = T)
    if (coords && !is.null(channels(x))) {
      df <- dplyr::left_join(df,
                             channels(x)[, c("electrode", "x", "y")],
                             by = "electrode")
    }
  }
  df
}



#' Check consistency of labels
#'
#' Internal function for checking 1) whether the labels submitted are a mixture
#' of hierarchical and non-hierarchical types 2) whether the labels submitted
#' are present in the data
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param cond_labs labels submitted by the user
#' @param data_labs labels from the actual data
#' @keywords internal

label_check <- function(cond_labs,
                        data_labs) {

  if (all(grepl("/", cond_labs))) {
    lab_check <- cond_labs %in% data_labs
    } else if (any(grepl("/",
                         cond_labs))) {
      stop("Do not mix hierarchical and non-hierarchical event labels.")
      } else {
        # Check if there is a hierarchical separator "/". If so,
        # split the labels
        if (any(grepl("/",
                      data_labs))) {
          split_labels <- strsplit(data_labs,
                                   "/")

          lab_check <- lapply(cond_labs,
                              function(x) vapply(split_labels,
                                                 function(i) x %in% i,
                                                 logical(1)))
          #condense to a single TRUE or FALSE for each label
          lab_check <- vapply(lab_check,
                              any,
                              logical(1))
        } else {
          lab_check <- cond_labs %in% data_labs
        }
      }
}

#' Convert to 3d matrix
#'
#' @param data data to be converted
#' @param ... additional parameters
#' @keywords internal
conv_to_mat <- function(data,...) {
  UseMethod("conv_to_mat", data)
}

conv_to_mat.default <- function(data, ...) {
  stop("Not implemented for objects of class", class(data))
}

#' @describeIn conv_to_mat Convert eeg_epochs to 3D matrix
conv_to_mat.eeg_epochs <- function(data, ...) {
  n_epochs <- length(unique(data$timings$epoch))
  n_channels <- ncol(data$signals)
  n_times <- length(unique(data$timings$time))
  data <- array(as.matrix(data$signals),
                dim = c(n_times, n_epochs, n_channels))
  data
}


is.true <- function(x) {
  !is.na(x) & x
}

#' Generate a zero centred vector
#'
#' @param vec_length how many points you want the vector to have
#' @importFrom stats median
#' @keywords internal
zero_vec <- function(vec_length) {
  out_seq <- seq(1, vec_length, by = 1)
  out_seq - stats::median(out_seq)
}

#' Pad a vector with zeros
#'
#' Pads a vector with zeros at the beginning and end.
#'
#' @param x vector to pad
#' @param n number of zeros to pad
#' @keywords internal
pad <- function(x, n) {c(rep(0, n), x, rep(0, n))}

#' Unpad a vector
#'
#' remove a specified number of points from the beginning and end of a vector
#' @param x Vector to unpad/trim
#' @param n Number of points to remove
#' @keywords internal
unpad <- function(x, n) {
  start <- n + 1
  end <- length(x) - n
  x <- x[start:end]
  x
}

#' Fix group delay
#'
#' Corrects a signal for the group delay of an FIR filter.
#'
#' @keywords internal
fix_grpdelay <- function(x, n, grp_delay) {
  start <- n + 1 + grp_delay
  end <- length(x) - n + grp_delay
  x <- x[start:end]
  x
}

#' Get number of samples
#' @param .data Object to get total number of sampling points from
#' @keywords internal
samples <- function(.data) {
  nrow(.data$signals)
}

#' Generate circle as radians
#' @keywords internal
circ_rad_fun <- function() {
  seq(0,
      2 * pi,
      length.out = 101)
}
