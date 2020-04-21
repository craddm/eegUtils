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


  if (length(dim(x$signals)) == 3) {
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
  } else {
    # Check which dimenstion to average over - participant or epoch
    if ("epoch" %in% x$dimensions) {
      avg_dim <- "epoch"
    } else if ("participant_id" %in% x$dimensions) {
      avg_dim <- "participant_id"
    }
    out_df <- aperm(x$signals, c("time",
                                 "frequency",
                                 avg_dim,
                                 "electrode"))
    dim_size <- dim(out_df)
    dim(out_df) <- c(dim_size[1] * dim_size[2] * dim_size[3],
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

  if ("participant_id" %in% x$dimensions) {
    out_df <- dplyr::left_join(out_df,
                               x$epochs,
                               by = c("participant_id", "epoch"))
   } else {
    # out_df$epoch <- as.numeric(out_df$epoch)
     out_df <- dplyr::left_join(out_df,
                                x$epochs,
                                by = c("epoch"))
   }
  if (long) {
    return(tidyr::gather(out_df,
                         electrode,
                         power,
                         channel_names(x)))
  }

  out_df
}



#' Convert `eeg_lm` to data.frame
#'
#' Convert an object of class \code{eeg_data} into a standard `data.frame`.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x Object of class \code{eeg_data}
#' @param row.names Kept for compatibility with S3 generic, ignored.
#' @param optional Kept for compatibility with S3 generic, ignored.
#' @param long Convert to long format. Defaults to FALSE.
#' @param coords Include electrode coordinates in output. Only possible when
#'   long = TRUE.
#' @param values Defaults to "coefficients", returning fitted coefficient values.
#' @param ... arguments for other as.data.frame commands
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
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

#    df$electrode <- as.character(df$electrode)

    if (coords && !is.null(channels(x))) {
      df <- dplyr::left_join(df,
                             channels(x)[, c("electrode", "x", "y")],
                             by = "electrode")
    }
  }
  df
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



