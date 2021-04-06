#' Function to create an S3 object of class "eeg_stats".
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param statistic Calculated statistic (e.g. t-statistic)
#' @param pvals calculated p-values for that statistic
#' @param chan_info String of character names for electrodes.
#' @param timings Unique timepoints remaining in the data.
#' @param method Type of statistical test
#' @keywords internal

eeg_stats <- function(statistic,
                      chan_info,
                      pvals,
                      timings,
                      method) {

  value <- list(statistic = statistic,
                pvals = tibble::as_tibble(pvals),
                chan_info = chan_info,
                timings = timings,
                method = method,
                version = utils::packageVersion("eegUtils"))
  class(value) <- "eeg_stats"
  value
}

new_eeg_lm <- function(coefficients,
                       std_err,
                       t_stats,
                       r_sq,
                       chan_info,
                       epochs,
                       timings,
                       formula) {

  stopifnot(is.data.frame(coefficients))
  stopifnot(is.data.frame(std_err))
  stopifnot(is.data.frame(t_stats))
  stopifnot(is.data.frame(r_sq))
  stopifnot(rlang::is_formula(formula))

  new_eeg_stats(
    coefficients = coefficients,
    std_err = std_err,
    t_stats = t_stats,
    r_sq = r_sq,
    timings = timings,
    epochs = epochs,
    chan_info = chan_info,
    formula = formula,
    class = "eeg_lm"
  )
}

new_eeg_stats <- function(...,
                          class = character()) {

  new_obj <- list(...)
  structure(
    new_obj,
    class = c(class, "eeg_stats")
  )
}

# new_eeg_tstats <- function(statistic) {
#   stopifnot(is.data.frame(statistic))
#   stopifnot(is.data.frame(timings))
#   #stopifnot(is.data.frame(t_stats))
#   #stopifnot(is.data.frame(r_sq))
#   #stopifnot(rlang::is_formula(formula))
#
#   new_eeg_stats(
#     statistic = statistic,
#     timings = timings,
#     epochs = epochs,
#     chan_info = chan_info,
#     formula = formula,
#     class = "eeg_tstats"
#   )
# }
