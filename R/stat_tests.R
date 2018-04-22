#' Compare epochs using an independent t-test
#'
#' @param data \code{eeg_epochs} object.
#' @param ... Other parameters passed to functions.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}

compare_epochs <- function(data, ...) {
  UseMethod("compare_epochs", data)
}

#' @importFrom purrr map_df
#' @param cond_label Conditions to compare
#' @param type Type of test to use. "1samp", "2samp"
#' @describeIn compare_epochs Compare differences across epochs

compare_epochs.eeg_epochs <- function(data, cond_label = NULL, type, ...) {

  if (identical(type, "1samp")) {
    data$signals <- split(data$signals, data$timings$time)
    data$signals <- purrr::map_df(data$signals, ~purrr::map(., calc_tstat))
  } else if (identical(type, "2samp")) {
  #  data$signals <- split(data$signals, list(data$timings$time, data$timings$)
  }
  eeg_stats(statistic = data$signals,
            chan_info = data$chan_info,
            timings = unique(data$timings$time))
}

#' Calculate one-sample t-statistic
#'
#' @param x vector of values to compare against mu
#' @param mu mean to compare x to.
#'

calc_tstat <- function(x, mu = 0) {
  tstat <- (mean(x) - mu) / (sd(x) / sqrt(length(x)))
}

#' Calculate two-sample independent t-statistic using Welch's formula
#'
#' @param x1 vector of values from one condition
#' @param x2 vector of values from one condition

calc_tstat_2 <- function(x1, x2) {
  tstat <- mean(x1) - mean(x2) / sqrt(var(x1) / length(x1) + var(x2) / length(x2))
}

#' Calculate p-value
#'
#' @param x t-stat
#' @param df Degrees of freedom
#' @param tails 1 or 2 for one- or two-tailed p-value.
#'

calc_pval <- function(x, df, tails = 2) {
  pval <- tails * pt(-abs(x), df)
}
