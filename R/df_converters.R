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

    df <- tidyr::pivot_longer(df,
                              cols = names(x$coefficients),
                              names_to = "electrode",
                              values_to = "amplitude")

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



