#' Compare epochs using an independent t-test
#'
#' @param data \code{eeg_epochs} object.
#' @param ... Other parameters passed to functions.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}

compare_epochs <- function(data, ...) {
  UseMethod("compare_epochs", data)
}

#' @importFrom tibble as_tibble
#' @param cond_label Conditions to test. A character vector of max length 2.
#' @param type Type of test to use. "1samp", "2samp"
#' @describeIn compare_epochs Compare differences across epochs

compare_epochs.eeg_epochs <- function(data, cond_label = NULL, type, ...) {

  if (identical(type, "1samp")) {
    #select only data specified in cond_label, if it exists
    if (!is.null(cond_label)) {
      data <- select_epochs(data, cond_label)
    }

    #matrix_t testing
    elecs <- colnames(data$signals)
    data$signals <- as.matrix(data$signals)
    n_epochs <- length(unique(data$timings$epoch))
    n_times <- length(unique(data$timings$time))
    dim(data$signals) <- c(n_times, n_epochs, ncol(data$signals))
    data$signals <- array_t(data$signals)
    colnames(data$signals) <- elecs
    data$signals <- tibble::as_tibble(data$signals)
    data$pvals <- apply(data$signals, 2,
                        function(x) calc_pval(x, df = (n_epochs - 1)))
  } else if (identical(type, "2samp")) {
    if (length(cond_label) != 2) {
      stop("Two condition labels must be supplied.")
    }

    elecs <- colnames(data$signals)

    # select epochs from data frame

    # warning - this does not validate that the epochs in each condition are
    # unique - it's possible for an epoch to contain more than one event, and
    # this selects based on *events*
    n_times <- length(unique(data$timings$time))
    cond1 <- select_epochs(data, cond_label[[1]])
    n_epochs <- length(unique(cond1$timings$epoch))
    deg_f <- n_epochs

    cond1 <- as.matrix(cond1$signals)
    dim(cond1) <- c(n_times, n_epochs, ncol(data$signals))

    cond2 <- select_epochs(data, cond_label[[2]])
    n_epochs <- length(unique(cond2$timings$epoch))
    deg_f <- deg_f + n_epochs

    cond2 <- as.matrix(cond2$signals)
    dim(cond2) <- c(n_times, n_epochs, ncol(data$signals))

    #cond1 <- split(cond1$signals, cond1$timings$time)
    #cond2 <- split(cond2$signals, cond2$timings$time)
    data$signals <- calc_tstat_2(cond1, cond2)
    colnames(data$signals) <- elecs
    data$signals <- tibble::as_tibble(data$signals)
    data$pvals <- apply(data$signals, 2,
                        function(x) calc_pval(x, df = (deg_f - 2)))
  }
  eeg_stats(statistic = data$signals,
            pvals = data$pvals,
            chan_info = data$chan_info,
            timings = unique(data$timings$time))
}

#' Array one-sample t-statistic
#'
#' @param x 3D array (timepoints x epochs x electrodes)
#' @param mu mean to test against
#' @importFrom matrixStats rowSds
#' @author Matt Craddock, \email{matt@@mattcraddock.com}

array_t <- function(x, mu = 0) {
  #calculate means for each combination of timepoint and electrode
  tp_means <- colMeans(aperm(x, c(2, 1, 3)))

  #calculate standard deviation for each combination of timepoint and electrode
  tp_sds <- vapply(1:dim(x)[[3]], function(i) matrixStats::rowSds(x[,,i]), numeric(dim(x)[[1]]))

  new_t <- (tp_means - mu) / (tp_sds /sqrt(dim(x)[[2]]))
  new_t
}


#' Calculate two-sample independent t-statistic using Welch's formula
#'
#' @param x1 vector of values from one condition
#' @param x2 vector of values from one condition
#' @author Matt Craddock, \email{matt@@mattcraddock.com}

calc_tstat_2 <- function(x1, x2) {


  x1_means <- colMeans(aperm(x1, c(2, 1, 3)))
  x2_means <- colMeans(aperm(x2, c(2, 1, 3)))

  x1_vars <- vapply(1:dim(x1)[[3]], function(i) matrixStats::rowVars(x1[, , i]), numeric(dim(x1)[[1]]))
  x2_vars <- vapply(1:dim(x2)[[3]], function(i) matrixStats::rowVars(x2[, , i]), numeric(dim(x2)[[1]]))

  #adapt to calculate df too...

  hmmz <- (x1_means - x2_means) / sqrt(x1_vars / dim(x1)[[2]] + x2_vars / dim(x2)[[2]])
  hmmz
}

#' Calculate p-value
#'
#' @param x t-stat
#' @param df Degrees of freedom
#' @param tails 1 or 2 for one- or two-tailed p-value.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}

calc_pval <- function(x, df, tails = 2) {
  pval <- tails * pt(-abs(x), df)
}
