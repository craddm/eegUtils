#' Compare epochs using an independent t-test
#'
#' @param data \code{eeg_epochs} object.
#' @param ... Other parameters passed to functions.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}

compare_epochs <- function(data, ...) {
  UseMethod("compare_epochs", data)
}

#' @describeIn compare_epochs Compare differences across epochs
compare_epochs.eeg_epochs <- function(data, ...) {

  data$signals <- split(data$signals, data$timings$time)
  data$signals <- purrr::map_df(data$signals, ~purrr::map(., calc_tstat))
  eeg_stats(statistic = data$signals,
            chan_info = data$chan_info,
            timings = unique(data$timings$time))
}

#' Calculate t-statistic
#'
#' @param x vector of values to compare against mu
#' @param mu mean to compare x to.
#'
calc_tstat <- function(x, mu = 0) {
  tstat <- (mean(x) - mu) / (sd(x) / sqrt(length(x)))
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
