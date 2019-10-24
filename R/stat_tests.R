#' Compare epochs using an independent t-test.
#'
#' @param data \code{eeg_epochs} object.
#' @param ... Other parameters passed to functions.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

compare_epochs <- function(data,
                           ...) {
  UseMethod("compare_epochs", data)
}

#' @importFrom tibble as_tibble
#' @param cond_label condition(s) to test.
#' @param type Type of test to use. "1samp", "2samp"
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @describeIn compare_epochs Compare differences across epochs
#' @keywords internal

compare_epochs.eeg_epochs <- function(data,
                                      cond_label = NULL,
                                      type, ...) {

  if (identical(type, "1samp")) {
    #select only data specified in cond_label, if it exists
    if (!is.null(cond_label)) {
      data <- select_epochs(data, cond_label)
    }

    #matrix_t testing
    elecs <- channel_names(data)
    data$signals <- as.matrix(data$signals)
    n_epochs <- length(unique(data$timings$epoch))
    n_times <- length(unique(data$timings$time))
    dim(data$signals) <- c(n_times, n_epochs, ncol(data$signals))
    data$signals <- array_t(data$signals)
    colnames(data$signals) <- elecs
    data$pvals <- calc_pval(data$signals,
                            df = (n_epochs - 1))
    data$signals <- tibble::as_tibble(data$signals)
    type <- "One-sample t-test"
  } else if (identical(type, "2samp")) {

    elecs <- channel_names(data)
    n_elecs <- length(elecs)
    n_times <- length(unique(data$timings$time))
    n_epochs <- length(unique(data$timings$epoch))
    deg_f <- n_epochs - 2

    data$signals <- as.data.frame(data)
    conds <- unique(data$signals[[cond_label]])

    if (length(conds) != 2) {
      stop("Condition must have two levels.")
    }

    cond1 <- data$signals[[cond_label]] == conds[1]
    cond2 <- data$signals[[cond_label]] == conds[2]

    #mat1 <- as.matrix(data$signals[cond1, elecs])
    #mat2 <- as.matrix(data$signals[cond2, elecs])
    #dim(mat1) <- c(n_times, sum(cond1) / n_times, n_elecs)
    #dim(mat2) <- c(n_times, sum(cond2) / n_times, n_elecs)

    data$signals <- calc_tstat_2(array(as.matrix(data$signals[cond1, elecs]),
                                       dim = c(n_times, sum(cond1) / n_times, n_elecs)),
                                 array(as.matrix(data$signals[!cond1, elecs]),
                                       dim = c(n_times, sum(cond1) / n_times, n_elecs)))
    colnames(data$signals) <- elecs
    data$pvals <- calc_pval(data$signals,
                            df = deg_f)
    data$signals <- tibble::as_tibble(data$signals)
    type <- "Welch two-sample t-test "
  }

  eeg_stats(statistic = data$signals,
            pvals = data$pvals,
            chan_info = data$chan_info,
            timings = unique(data$timings$time),
            method = type)
}

#' Array one-sample t-statistic
#'
#' @param x 3D array (timepoints x epochs x electrodes)
#' @param mu mean to test against
#' @importFrom matrixStats rowSds
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

array_t <- function(x, mu = 0) {
  #calculate means for each combination of timepoint and electrode
  tp_means <- colMeans(aperm(x,
                             c(2, 1, 3)))

  #calculate standard deviation for each combination of timepoint and electrode
  tp_sds <- vapply(1:dim(x)[[3]],
                   function(i) matrixStats::rowSds(x[, , i]),
                   numeric(dim(x)[[1]]))

  new_t <- (tp_means - mu) / (tp_sds / sqrt(dim(x)[[2]]))
  new_t
}

#' Calculate two-sample independent t-statistic using Welch's formula
#'
#' @param x1 matrix of values from one condition
#' @param x2 matrix of values from one condition
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

calc_tstat_2 <- function(x1, x2) {

  x1_means <- colMeans(aperm(x1, c(2, 1, 3)))
  x2_means <- colMeans(aperm(x2, c(2, 1, 3)))

  x1_vars <- vapply(1:dim(x1)[[3]],
                    function(i) matrixStats::rowVars(x1[, , i]),
                    numeric(dim(x1)[[1]]))
  x2_vars <- vapply(1:dim(x2)[[3]],
                    function(i) matrixStats::rowVars(x2[, , i]),
                    numeric(dim(x2)[[1]]))

  #adapt to calculate df too...
  #welch formula
  t_stats <- (x1_means - x2_means) /
    sqrt(x1_vars / dim(x1)[[2]] +
           x2_vars / dim(x2)[[2]])
  t_stats
}

#' Calculate p-value
#'
#' @param x Vector/matrix of t-stats
#' @param df Degrees of freedom.
#' @param tails 1 or 2 for one- or two-tailed p-value.
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

calc_pval <- function(x, df, tails = 2) {
  pval <- tails * pt(-abs(x), df)
}

#' Calculate F-statistic for independent samples
#'
#' @param x1 condition one
#' @param x2 condition two
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @keywords internal

calc_fstat <- function(x1, x2) {
  x1_means <- colMeans(aperm(x1,
                             c(2, 1, 3)))
  x2_means <- colMeans(aperm(x2,
                             c(2, 1, 3)))

  grand_means <- (x1_means + x2_means) / 2

  x1_test <- vapply(1:dim(x1)[[3]],
                    function(i) x1[, , i] - grand_means[, i],
                    matrix(0, nrow = nrow(x1),
                           ncol = ncol(x1)))

  x2_test <- vapply(1:dim(x2)[[3]],
                    function(i) x2[, , i] - grand_means[, i],
                    matrix(0, nrow = nrow(x2),
                           ncol = ncol(x2)))

  x3_test_sst <- vapply(1:dim(x1)[[3]],
                        function(i) cbind(x1_test[, , i],
                                          x2_test[, , i]) ^ 2,
                        matrix(1,
                               nrow = nrow(x1_test),
                               ncol = ncol(x1_test) + ncol(x2_test)))

  x3_test_sst <- apply(x3_test_sst,
                       c(1, 3),
                       sum)
  x3_ssm <- ncol(x1) * (x1_means - grand_means) ^ 2 +
    ncol(x2) * (x2_means - grand_means) ^ 2

  x3_ssr <- x3_test_sst - x3_ssm

  x4_msm <- x3_ssm / 1
  x4_msr <- x3_ssr / (ncol(x1) + ncol(x2) - 2)
  x4_f <- x4_msm / x4_msr
  x4_f
}



permute_two <- function(data,
                        cond1,
                        n_times,
                        n_elecs) {
  perms <- sample(cond1)
  calc_tstat_2(array(data[perms, ],
                     dim = c(n_times, sum(perms)/n_times, n_elecs)),
               array(data[!perms, ],
                     dim = c(n_times, sum(perms)/n_times, n_elecs)))
}
