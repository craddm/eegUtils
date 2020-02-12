#' Simple absolute value thresholding
#'
#' Reject data based on a simple absolute amplitude threshold. This marks any timepoint
#' from any electrode.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data} or \code{eeg_epochs}.
#' @param threshold In microvolts. If one value is supplied, it will be treated
#'   as a +- value.
#' @param reject If TRUE, remove marked data immediately, otherwise mark for
#'   inspection/rejection. Defaults to FALSE.
#' @examples
#' ar_thresh(demo_epochs, c(100))
#' @export

ar_thresh <- function(data,
                      threshold,
                      reject = FALSE) {
  UseMethod("ar_thresh", data)
}

#' @describeIn ar_thresh Reject data using a simple threshold.
#' @export
ar_thresh.eeg_data <- function(data,
                               threshold,
                               reject = FALSE) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- check_thresh(data, threshold)
  crossed_thresh <- rowSums(crossed_thresh) == 0

  if (reject) {
    message("Removing ",
            sum(!crossed_thresh), " (",
            round(sum(!crossed_thresh) / nrow(data$signals) * 100, 2),
            "%) timepoints.")
    data$reject$timings <- data$timings[crossed_thresh, ]
    data$timings <- data$timings[crossed_thresh, ]
    data$events <- data$events[data$events$event_time %in% data$timings$time, ]

    data$signals <- data$signals[crossed_thresh, ]

  } else {
    data$reject$timings <- data$timings[crossed_thresh, ]
  }
  data
}

#' @describeIn ar_thresh Reject data using a simple threshold.
#' @export
ar_thresh.eeg_epochs <- function(data,
                                 threshold,
                                 reject = FALSE) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- check_thresh(data, threshold)

  crossed_thresh <- rowSums(crossed_thresh) == 1
  rej_epochs <- unique(data$timings$epoch[crossed_thresh])
  message(paste(length(rej_epochs), "epochs contain samples above threshold."))
  if (reject) {
    message("Removing ", length(rej_epochs), " epochs.")
    data <- select_epochs(data,
                          epoch_no = rej_epochs,
                          keep = FALSE)
  } else {
    data$reject$epochs <- rej_epochs
    data$reject$timings <- data$timings[crossed_thresh,]
  }
  data
}

check_thresh <- function(data, threshold) {
  crossed_thresh <- data$signals > max(threshold) |
    data$signals < min(threshold)

  upper_thresh <- data$signals > max(threshold)
  lower_thresh <- data$signals < min(threshold)
  total_data <- prod(dim(data$signals))

  message(sum(upper_thresh),
          " (", round(sum(upper_thresh)/total_data * 100, 2), "%) ",
          "samples above ", max(threshold) , " uV threshold.")
  message(sum(lower_thresh),
          " (", round(sum(lower_thresh)/total_data * 100, 2), "%) ",
          "samples below ", min(threshold) , " uV threshold.")
  crossed_thresh
}

#' Channel statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param ... Other parameters passed to the functions.
#' @keywords internal

channel_stats <- function(data, ...) {
  UseMethod("channel_stats", data)
}

#' @describeIn channel_stats Calculate channel statistics for \code{eeg_data}
#'   objects.
#' @keywords internal
channel_stats.eeg_data <- function(data, ...) {

  chan_means <- colMeans(data$signals)
  chan_sds <- apply(data$signals, 2, stats::sd)
  chan_var <- chan_sds^2
  #chan_kurt <- as.numeric(as.data.table(data$signals)[, lapply(.SD, kurtosis)])
  chan_kurt <- apply(data$signals, 2, kurtosis)

  data.frame(electrode = names(data$signals),
             means = chan_means,
             sds = chan_sds,
             variance = chan_var,
             kurtosis = chan_kurt
             )
}

#' Epoch statistics
#'
#' Calculate various statistics for each epoch in the data
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param ... Other parameters passed to the functions.
#' @keywords internal

epoch_stats <- function(data, ...) {
  UseMethod("epoch_stats", data)
}

#' @describeIn epoch_stats Calculate statistics for each epoch.
#' @keywords internal
epoch_stats.eeg_epochs <- function(data, ...) {
  data$signals$epoch <- data$timings$epoch
  data <- data.table::data.table(as.data.frame(data$signals))
  epoch_vars <- data[, lapply(.SD, var), by = epoch]
  epoch_kur <- data[, lapply(.SD, kurtosis), by = epoch]
  epoch_max <- data[, lapply(.SD, max), by = epoch]
  epoch_min <- data[, lapply(.SD, min), by = epoch]
  stats_out <- data.table::rbindlist(list(max = epoch_max,
                                          min = epoch_min,
                                          variance = epoch_vars,
                                          kurtosis = epoch_kur),
                                     idcol = "measure")
  stats_out
}

#' Calculate kurtosis
#'
#' @param data Data to calculate kurtosis for
#' @keywords internal

kurtosis <- function(data) {
  m4 <- mean((data - mean(data)) ^ 4)
  kurt <- m4 / (stats::sd(data) ^ 4) - 3
  kurt
}

#' Remove EOG using regression
#'
#' Calculates and removes the contribution of eye movements to the EEG signal using
#' least-squares regression.
#'
#' @param data data to regress - \code{eeg_data} or \code{eeg_epochs}
#' @param heog Horizontal EOG channel labels
#' @param veog Vertical EOG channel labels
#' @param bipolarize Bipolarize the EOG channels. Only works when four channels
#'   are supplied (2 HEOG and 2 VEOG).
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @export

ar_eogreg <- function(data,
                      heog,
                      veog,
                      bipolarize = TRUE) {
  UseMethod("ar_eogreg", data)
}

#' @rdname ar_eogreg
#' @export
ar_eogreg.eeg_data <- function(data,
                               heog,
                               veog,
                               bipolarize = TRUE) {

  eogreg(data,
         heog,
         veog,
         bipolarize)
}

#' @rdname ar_eogreg
#' @export
ar_eogreg.eeg_epochs <- function(data,
                                 heog,
                                 veog,
                                 bipolarize = TRUE) {

  eogreg(data,
         heog,
         veog,
         bipolarize)

}

#' @noRd
eogreg <- function(data,
                   heog,
                   veog,
                   bipolarize) {

  if (bipolarize) {
    EOG <- bip_EOG(data$signals, heog, veog)
  } else {
    HEOG <- data$signals[, heog, drop = TRUE]
    VEOG <- data$signals[, veog, drop = TRUE]
    EOG <- data.frame(HEOG, VEOG)
  }

  data_chans <- channel_names(data)[!channel_names(data) %in% c(heog, veog)]
  hmz <- solve(crossprod(as.matrix(EOG)),
               crossprod(as.matrix(EOG),
                         as.matrix(data$signals[, data_chans])))
  data$signals[, data_chans] <- data$signals[, data_chans] - crossprod(t(as.matrix(EOG)), hmz)
  data
}

#' @noRd
bip_EOG <- function(data,
                    HEOG,
                    VEOG) {
  HEOG <- data[, HEOG[1]] - data[, HEOG[2]]
  VEOG <- data[, VEOG[1]] - data[, VEOG[2]]
  EOG <- data.frame(HEOG, VEOG)
  EOG
}

#' Detect high correlation with EOG channels
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param decomp ICA decomposition
#' @param ... Other parameters
#' @export

ar_eogcor <- function(decomp, ...) {
  UseMethod("ar_eogcor", decomp)
}

#' @param data Original data
#' @param HEOG Horizontal eye channels
#' @param VEOG Vertical eye channels
#' @param threshold Threshold for correlation (r). Defaults to NULL,
#'   automatically determining a threshold.
#' @param plot Plot correlation coefficient for all components
#' @describeIn ar_eogcor Method for eeg_ICA objects.
#' @export
ar_eogcor.eeg_ICA <- function(decomp,
                              data,
                              HEOG,
                              VEOG,
                              threshold = NULL,
                              plot = TRUE,
                              ...) {

  if (!is.null(threshold)) {
    if (threshold > 1 | threshold < 0) {
      stop("Threshold must be between 0 and 1.")
    }
    heog_threshold <- threshold
    veog_threshold <- threshold
  }


  EOG_corrs <- abs(stats::cor(decomp$signals,
                              bip_EOG(data$signals,
                                      HEOG,
                                      VEOG)))


  if (is.null(threshold)) {
    mean_corrs <- colMeans(EOG_corrs)
    sd_corrs <- matrixStats::colSds(EOG_corrs)
    heog_threshold <- mean_corrs[[1]] + 4 * sd_corrs[[1]]
    veog_threshold <- mean_corrs[[2]] + 4 * sd_corrs[[2]]
    message("Estimated HEOG threshold: ", round(heog_threshold, 2))
    message("Estimated VEOG threshold: ", round(veog_threshold, 2))
  }

  if (plot) {
    graphics::par(mfrow = c(1, 2))
    plot(EOG_corrs[, 1])
    graphics::abline(h = heog_threshold)
    plot(EOG_corrs[, 2])
    graphics::abline(h = veog_threshold)
  }

  crossed_thresh <- cbind(EOG_corrs[, 1] > heog_threshold,
                          EOG_corrs[, 2] > veog_threshold)
  above_thresh <- apply(crossed_thresh,
                        1, any)
  above_thresh <- channel_names(decomp)[above_thresh]
  above_thresh
}

#' Detect low autocorrelation of ICA components
#'
#' Low autocorrelation can be a sign of a poor quality channel or component.
#' Often these are noisy, poor contact, or heavily contaminated with muscle
#' noise. Low autocorrelation at a lag of 20ms is often associated with muscle
#' noise.
#'
#' @param data \code{eeg_ICA} object
#' @param ... additional parameters
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @references Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide
#'   to the selection of independent components of the electroencephalogram for
#'   artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
#'   10.1016/j.jneumeth.2015.02.025
#' @export
ar_acf <- function(data, ...) {
  UseMethod("ar_acf", data)
}

#' @param ms Time lag to check ACF, in milliseconds. Defaults to 20 ms.
#' @param plot Produce plot showing ACF and threshold for all EEG components.
#' @param verbose Print informative messages. Defaults to TRUE.
#' @param threshold Specify a threshold for low ACF. NULL estimates the threshold automatically.
#' @describeIn ar_acf Autocorrelation checker for \code{eeg_ICA} objects
#' @importFrom stats sd cor
#' @importFrom graphics abline par
#' @export
ar_acf.eeg_ICA <- function(data,
                           ms = 20,
                           plot = TRUE,
                           verbose = TRUE,
                           threshold = NULL,
                           ...) {

  time_lag <- round(data$srate * (ms/1000))
  low_acf <- apply(data$signals, 2,
                   function(x) stats::acf(x, time_lag, plot = FALSE)$acf[time_lag + 1, 1, 1])
  if (is.null(threshold)) {
    threshold <- mean(low_acf) - 2 * stats::sd(low_acf)
  }
  if (verbose) {
    message("Estimating autocorrelation at ", ms, "ms lag.")
    message("Estimated ACF threshold: ", round(threshold, 2))
  }
  if (plot) {
    plot(low_acf)
    graphics::abline(h = threshold)
  }
  low_acf <- channel_names(data)[low_acf < threshold]
  message("Subthreshold components: ", paste0(low_acf, sep = " "))
  low_acf
}


#' Detect high channel focality of ICA components
#'
#' Detect components that load heavily on a single channel. Looks for components
#' that have one particular channel that has a particularly high z-score.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data \code{eeg_ICA} object
#' @param plot Produce plot showing max z-scores and threshold for all ICA
#'   components.
#' @param threshold Specify a threshold for high focality. NULL estimates the
#'   threshold automatically.
#' @param verbose Print informative messages.
#' @param measure Use maximum "max" or "kurtosis".
#' @param ...  additional parameters
#' @references Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide
#'   to the selection of independent components of the electroencephalogram for
#'   artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
#'   10.1016/j.jneumeth.2015.02.025
#' @export
ar_chanfoc <- function(data,
                       plot = TRUE,
                       threshold = NULL,
                       verbose = TRUE,
                       measure = "max",
                       ...) {

  n_comps <- length(channel_names(data))
  zmat <- apply(data$mixing_matrix[, 1:n_comps], 2, scale)
  if (identical(measure, "max")) {
    max_weights <- apply(zmat, 2, function(x) max(abs(x)))
  } else if (identical(measure, "kurtosis")) {
    max_weights <- apply(zmat, 2, kurtosis)
  }
  if (is.null(threshold)) {
    threshold <- mean(max_weights) + 2 * stats::sd(max_weights)
    if (verbose) {
      message("Estimated threshold: ", round(threshold, 2))
    }
  }

  if (plot) {
    graphics::plot(max_weights)
    graphics::abline(h = threshold)
  }

  chan_foc <- channel_names(data)[max_weights > threshold]

  if (verbose) {
    message("Components with high channel focality: ", paste0(chan_foc, sep = " "))
  }

  chan_foc
}

#' Detect high trial focality of ICA components
#'
#' Detect components that load heavily on a small number of trials. Looks for components
#' that have one particular trial that has a particularly high z-score.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data \code{eeg_ICA} object
#' @param plot Produce plot showing max z-scores and threshold for all ICA
#'   components.
#' @param threshold Specify a threshold (z-score) for high focality. NULL estimates the
#'   threshold automatically.
#' @param verbose Print informative messages.
#' @references Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide
#'   to the selection of independent components of the electroencephalogram for
#'   artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
#'   10.1016/j.jneumeth.2015.02.025
#' @export
ar_trialfoc <- function(data,
                        plot = TRUE,
                        threshold = NULL,
                        verbose = TRUE) {
  n_comps <- length(channel_names(data))
  n_epochs <- nrow(epochs(data))
  n_times <- nrow(data$signals) / n_epochs
  zmat <- as.matrix(data$signals)
  dim(zmat) <- c(n_times, n_epochs, n_comps)
  zmat <- apply(zmat, c(2, 3), function(x) diff(range(x)))
  zmat <- abs(apply(zmat, 2, scale))

  if (is.null(threshold)) {
    threshold <- mean(matrixStats::colMaxs(zmat)) + 2 * stats::sd(matrixStats::colMaxs(zmat))
    if (verbose) {
      message("Estimated trial focality threshold (z): ",
              round(threshold, 2))
    }
  }

  if (plot) {
    plot(matrixStats::colMaxs(zmat))
    graphics::abline(h = threshold)
  }

  channel_names(data)[matrixStats::colMaxs(zmat) > threshold]
}

