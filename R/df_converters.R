#' Convert `eeg_data` to `data.frame`
#'
#' Convert an object of class `eeg_data` into a standard data.frame.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_data`
#' @param row.names Kept for compatibility with S3 generic, ignored.
#' @param optional Kept for compatibility with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param events Include events in output.
#' @param coords Include electrode coordinates in output. Only possible when
#'   long = TRUE.
#' @param ... arguments for other as.data.frame commands
#' @export

as.data.frame.eeg_data <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   long = FALSE,
                                   events = FALSE,
                                   coords = FALSE,
                                   ...) {
  chan_names <- channel_names(x)
  new_df <- data.frame(x$signals,
                       x$timings)

  if (long) {
    new_df <- tidyr::pivot_longer(
      new_df,
      cols = tidyr::all_of(chan_names),
      names_to = "electrode",
      values_to = "amplitude"
      )
  }

  if (events) {
    new_df <- dplyr::left_join(new_df,
                           x$events,
                           by = c("sample" = "event_onset"))
  }
  new_df
}


#' Convert `eeg_tfr` objects to a `data.frame`
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_tfr`
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param ... arguments for other `as.data.frame` commands
#' @return A `data.frame` or `tibble`
#' @export
as.data.frame.eeg_tfr <- function(x,
                                  row.names = NULL,
                                  optional = FALSE,
                                  long = FALSE,
                                  ...) {

  if (length(dim(x$signals)) == 3) {
    # Note, as of 0.5.0.9000, eeg_tfr always has more than 3 dimensions, so this
    # is for compatibility
    out_df <- aperm(x$signals,
                    c("time",
                      "frequency",
                      "electrode"))
    dim_size <- dim(out_df)
    dim(out_df) <- c(dim_size[1] * dim_size[2],
                     dim_size[3])
    out_df <- as.data.frame(out_df)
    names(out_df) <- dimnames(x$signals)$electrode
    out_df$time <- rep(as.numeric(dimnames(x$signals)$time),
                       dim_size[2])
    out_df$frequency <- rep(as.numeric(dimnames(x$signals)$frequency),
                            each = dim_size[1])
    out_df$epoch <- unique(x$epochs$epoch)
  } else if (inherits(x, "eeg_group")) {
      out_df <- aperm(x$signals,
                      c("time",
                        "frequency",
                        "epoch",
                        "participant_id",
                        "electrode"))
      dim_size <- dim(out_df)
      dim(out_df) <- c(prod(dim_size[1:4]),
                       dim_size[5])
      colnames(out_df) <- dimnames(x$signals)$electrode
      out_df <- as.data.frame(out_df)
      out_df$time <- rep(as.numeric(dimnames(x$signals)$time),
                         prod(dim_size[2:4]))
      out_df$frequency <- rep(as.numeric(dimnames(x$signals)$frequency),
                              each = dim_size[1])
      out_df$epoch <- rep(x$epochs$epoch, each =
                            dim_size[1] * dim_size[2])
      out_df$participant_id <- rep(x$epochs$participant_id,
                                   each = dim_size[1] * dim_size[2])

    } else {
      out_df <- aperm(x$signals,
                      c("time",
                        "frequency",
                        "epoch",
                        "electrode"))
      dim_size <- dim(out_df)
      dim(out_df) <- c(prod(dim_size[1:3]),
                       dim_size[4])
      out_df <- as.data.frame(out_df)
      names(out_df) <- dimnames(x$signals)$electrode
      out_df$time <- rep(as.numeric(dimnames(x$signals)$time),
                         dim_size[3] * dim_size[2])#out_df$time)
      out_df$frequency <- rep(as.numeric(dimnames(x$signals)$frequency),
                              each = dim_size[1])
      out_df$epoch <- rep(x$epochs$epoch, each =
                            dim_size[1] * dim_size[2])
      out_df$participant_id <- rep(x$epochs$participant_id,
                                   each = dim_size[1] * dim_size[2])
      }

  if ("participant_id" %in% colnames(out_df)) {
    out_df <- dplyr::left_join(out_df,
                               x$epochs,
                               by = c("participant_id", "epoch"))
  } else {
    out_df <- dplyr::left_join(out_df,
                               x$epochs,
                               by = "epoch")
  }
  if (long) {
    return(tidyr::pivot_longer(out_df,
                               cols = channel_names(x),
                               names_to = "electrode",
                               values_to = "power"))
  }
  out_df
}



#' Convert `eeg_lm` to `data.frame`
#'
#' Convert an object of class `eeg_data` into a standard `data.frame`.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_lm`
#' @param row.names Kept for compatibility with S3 generic, ignored.
#' @param optional Kept for compatibility with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param coords Include electrode coordinates in output. Only possible when
#'   long = TRUE.
#' @param values Defaults to "coefficients", returning fitted coefficient values.
#' @param ... arguments for other as.data.frame commands
#' @export
as.data.frame.eeg_lm <- function(x,
                                 row.names = NULL,
                                 optional = FALSE,
                                 long = FALSE,
                                 coords = TRUE,
                                 values = c("coefficients",
                                            "std_err",
                                            "t_stats",
                                            "r_sq"),
                                 ...) {

  values <- match.arg(values)

  if (identical(values, "r_sq")) {
    df <- data.frame(x$r_sq)
  } else {
    df <- data.frame(x[[values]],
                     time = x$timings$time,
                     epoch = x$timings$epoch,
                     stringsAsFactors = FALSE)
    if (!is.null(x$epochs)) {
      df <- dplyr::left_join(df,
                             x$epochs,
                             by = "epoch")
    }
  }

  if (long) {

    val_name <- switch(values,
                       "coefficients" = "amplitude",
                       "t_stats" = "statistic",
                       "std_err" = "std_err",
                       "r_sq" = "r_sq")

    df <- tidyr::pivot_longer(df,
                              cols = names(x$coefficients),
                              names_to = "electrode",
                              values_to = val_name)

    if (coords && !is.null(channels(x))) {
      df <- dplyr::left_join(df,
                             channels(x)[, c("electrode", "x", "y")],
                             by = "electrode")
    }
  }
  df
}

#' Convert `eeg_stats` objects to data frames
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_stats`
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param coords Include electrode coordinates in output (ignored if long = FALSE)
#' @param ... arguments for other as.data.frame commands
#'
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
                        factor_key = TRUE)
    if (coords && !is.null(channels(x))) {
      df <- dplyr::left_join(df,
                             channels(x)[, c("electrode", "x", "y")],
                             by = "electrode")
    }
  }
  df
}

#' Convert `eeg_epochs` object to data.frame
#'
#' Convert an `eeg_epochs` object to a data.frame for use with whatever
#' packages you desire.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_epochs`
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

  chan_names <- channel_names(x)
  df <- data.frame(x$signals,
                   time = x$timings$time,
                   epoch = x$timings$epoch,
                   stringsAsFactors = FALSE)

  names(df)[1:length(chan_names)] <- chan_names

  if (long) {
    # When converting to long, use factor_key to keep channels in same order,
    # then convert back to character.
    df <- tidyr::pivot_longer(df,
                              names_to = "electrode",
                              values_to = "amplitude",
                              cols = channel_names(x))

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

#' Convert `eeg_evoked` object to data frame
#
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_evoked`
#' @param row.names Kept for compatability with S3 generic, ignored.
#' @param optional Kept for compatability with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param coords include electrode coordinates in output. Ignored if long = FALSE.
#' @param ... arguments for other as.data.frame commands
#'
#' @export

as.data.frame.eeg_evoked <- function(x,
                                     row.names = NULL,
                                     optional = FALSE,
                                     long = FALSE,
                                     coords = TRUE,
                                     ...) {

  is_group_df <- inherits(x, "eeg_group")
  chan_info <- channels(x)
  chan_names <- channel_names(x)
  epoch_info <- epochs(x)
  x <- cbind(x$signals,
             x$timings)

  if (is_group_df) {
    unique_cols <-
      c(names(epoch_info)[!(names(epoch_info) %in% names(x))],
        "participant_id", "epoch")
    x <- dplyr::left_join(x,
                          epoch_info[unique_cols],
                          by = c("participant_id",
                                 "epoch"))
  } else {
    x <- dplyr::left_join(x,
                          epoch_info,
                          by = "epoch")
  }

  if (long) {
    x <- tidyr::pivot_longer(x,
                             cols = dplyr::all_of(chan_names),
                             names_to = "electrode",
                             values_to = "amplitude")

    if (coords && !is.null(chan_info)) {
       x <- dplyr::left_join(x,
                             chan_info[, c("electrode", "x", "y")],
                             by = "electrode")
    }
  }

  x
}

#' Convert `eeg_ICA` object to data frame
#
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class `eeg_ICA`
#' @param row.names Kept for compatibility with S3 generic, ignored.
#' @param optional Kept for compatibility with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE
#' @param cond_labels add condition labels to data frame. Deprecated.
#' @param mixing If TRUE, outputs the mixing matrix. If FALSE, outputs source activations.
#' @param coords Adds electrode coordinates if TRUE; only if long data and the mixing matrix are requested.
#' @param ... arguments for other as.data.frame commands
#'
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
    df <- tidyr::pivot_longer(df,
                              names_to = "component",
                              values_to = "amplitude",
                              cols = dplyr::all_of(channel_names(x)))
  }
  df
}


