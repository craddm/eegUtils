#' Convert eeg_data to data.frame
#'
#' Convert an object of class \code{eeg_data} into a standard data.frame / tibble
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_data}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param events Include events in output.
#' @param ... arguments for other as.data.frame commands
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_data <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   long = FALSE,
                                   events = FALSE,
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
#' @param events Include events in output. Defaults to FALSE.
#' @param conditions Names of columns from the events table to include in output
#' @param cond_labels Add column tagging epochs with events that have matching
#'   labels.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_epochs <- function(x, row.names = NULL,
                                     optional = FALSE,
                                     long = FALSE,
                                     events = FALSE,
                                     conditions = NULL,
                                     cond_labels = NULL, ...) {

  if (!is.null(conditions)) {
    cond_cols <- names(events(x)) %in% c("epoch", conditions)
    if (!any(cond_cols)) {
      message("Columns not found.")
    }

    missing_rows <- apply(events(x)[, conditions],
                          2,
                          is.na)
    missing_rows <- apply(missing_rows, 1, any)

    events(x) <- events(x)[!missing_rows, cond_cols]
  }

  if (!is.null(cond_labels)) {
    lab_check <- label_check(cond_labels,
                             unique(list_epochs(x))$event_label)
    if (!all(lab_check)) {
      stop("Not all labels found. Use list_events to check labels.")
    }

    df <- lapply(seq_along(cond_labels),
                 function(ix) {
                   out_df <- as.data.frame(select_epochs(x,
                                                         cond_labels[[ix]]))
                   out_df$conditions <- cond_labels[[ix]]
                   out_df})

    df <- do.call("rbind", df)

    if (long) {
      df <- tidyr::gather(df,
                          electrode,
                          amplitude,
                          -time,
                          -epoch,
                          -conditions,
                          factor_key = TRUE)
    }

  } else {
    df <- data.frame(x$signals,
                     time = x$timings$time,
                     epoch = x$timings$epoch)
    # combine the new data frame with any condition labels from the events table
    if ("event_label" %in% names(x$events)) {
      df <- merge(df,
                  x$events[c("epoch", "event_label")],
                  by = "epoch")
      names(df)[which(names(df) == "event_label")] <- "conditions"
    }

    if (long) {
      df <- tidyr::gather(df,
                          electrode,
                          amplitude,
                          names(x$signals),
                          factor_key = TRUE)
    }
  }

  if (events) {
    df <- dplyr::left_join(df,
                           x$events,
                           by = c("sample" = "event_onset"))
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
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_evoked <- function(x, row.names = NULL,
                                     optional = FALSE, long = FALSE, ...) {
  if (class(x$signals) == "list") {
    cond_labels <- names(x$signals)
    df <- lapply(seq_along(cond_labels), function(ix) {
      out_df <- data.frame(x$signals[[ix]],
                           time = x$timings$time,
                           conditions = cond_labels[[ix]])
      out_df
    })
    df <- do.call("rbind", df)
    if (long) {
      df <- tidyr::gather(df,
                          "electrode",
                          "amplitude",
                          -time,
                          -conditions,
                          factor_key = T)
    }
  } else {
    df <- data.frame(x$signals, time = x$timings$time)
    if (long) {
      df <- tidyr::gather(df,
                          "electrode",
                          "amplitude",
                          -time,
                          factor_key = T)
    }
  }
  df
}

#' Convert \code{eeg_ICA} object to data frame
#
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @param x Object of class \code{eeg_ICA}
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param cond_labels add condition labels to data frame.
#' @param ... arguments for other as.data.frame commands
#'
#' @importFrom tidyr gather
#' @export

as.data.frame.eeg_ICA <- function(x,
                                  row.names = NULL,
                                  optional = FALSE,
                                  long = FALSE,
                                  cond_labels = NULL,
                                  ...) {

  if (is.null(x$comp_activations)) {
    x$comp_activations <- x$signals
  } else {
    x$comp_activations <- as.data.frame(x$comp_activations)
  }

  if (!is.null(cond_labels)) {
    lab_check <- label_check(cond_labels,
                             unique(list_epochs(x))$event_label)
    if (!all(lab_check)) {
      stop("Not all labels found. Use list_events to check labels.")
    }

    df <- lapply(seq_along(cond_labels), function(ix) {
      out_df <- as.data.frame(select_epochs(x, cond_labels[[ix]]))
      out_df$conditions <- cond_labels[[ix]]
      out_df
    })

    df <- do.call("rbind", df)

    if (long) {
      df <- tidyr::gather(df,
                          electrode,
                          amplitude,
                          -time,
                          -epoch,
                          -conditions,
                          factor_key = T)
    }
  } else {
    df <- data.frame(x$comp_activations, x$timings)
    df$sample <- NULL
    if (long) {
      df <- tidyr::gather(df,
                          electrode,
                          amplitude,
                          -time,
                          -epoch,
                          factor_key = T)
    }
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

  out_df <- as.data.frame.table(x$signals,
                                stringsAsFactors = FALSE)
  names(out_df) <- c(x$dimensions, "power")
  out_df$time <- as.numeric(out_df$time)
  out_df$frequency <- as.numeric(out_df$frequency)
  if (!long) {
    out_df <- tidyr::spread(out_df,
                            electrode,
                            power)
    return(out_df)
  }
  out_df
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

pad <- function(x, n) {x <- c(rep(0, n), x, rep(0, n))}

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
