#' Detect high component correlation with eye channels
#'
#' Checks the correlation between individual components of an `eeg_ICA`
#' decomposition and the electrooculogram channels of an `eeg_epochs` dataset.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param decomp ICA decomposition
#' @param data Original data
#' @param ... Other parameters
#' @references Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide
#'   to the selection of independent components of the electroencephalogram for
#'   artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
#'   10.1016/j.jneumeth.2015.02.025
#' @return A character vector of component names that break the threshold.
#' @export

ar_eogcor <- function(decomp,
                      data,
                      ...) {
  UseMethod("ar_eogcor", decomp)
}


#' @param HEOG Horizontal eye channels
#' @param VEOG Vertical eye channels
#' @param threshold Threshold for correlation (r). Defaults to NULL,
#'   automatically determining a threshold.
#' @param plot Plot correlation coefficient for all components
#' @param bipolarize Bipolarize the HEOG and VEOG channels?
#' @param verbose Print informative messages. Defaults to TRUE.
#' @describeIn ar_eogcor Method for eeg_ICA objects.
#' @export
ar_eogcor.eeg_ICA <- function(decomp,
                              data,
                              HEOG,
                              VEOG,
                              threshold = NULL,
                              plot = TRUE,
                              bipolarize = TRUE,
                              verbose = TRUE,
                              ...) {

  if (!is.null(threshold)) {
    if (threshold > 1 || threshold < 0) {
      stop("Threshold must be between 0 and 1.")
    }
    heog_threshold <- threshold
    veog_threshold <- threshold
  }

  EOG_corrs <- abs(stats::cor(decomp$signals,
                              if (bipolarize) {
                                bip_EOG(data$signals,
                                        HEOG,
                                        VEOG)
                              } else {
                                data$signals[, c("HEOG", "VEOG")]
                              }))

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
  if (verbose) {
    message("Components with high EOG correlation: ", paste0(above_thresh, sep = " "))
  }
  above_thresh
}

#' Detect low autocorrelation of ICA components
#'
#' Low autocorrelation can be a sign of a poor quality channel or component.
#' Often these are noisy, poor contact, or heavily contaminated with muscle
#' noise. Low autocorrelation at a lag of 20ms is often associated with muscle
#' noise.
#'
#' @param data `eeg_ICA` object
#' @param ... additional parameters
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @references Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide
#'   to the selection of independent components of the electroencephalogram for
#'   artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
#'   10.1016/j.jneumeth.2015.02.025
#' @examples
#' demo_sobi <- run_ICA(demo_epochs, pca = 10)
#' ar_acf(demo_sobi)
#' @return A character vector of component names that break the threshold.
#' @export
ar_acf <- function(data, ...) {
  UseMethod("ar_acf", data)
}

#' @param ms Time lag to check ACF, in milliseconds. Defaults to 20 ms.
#' @param plot Produce plot showing ACF and threshold for all EEG components.
#' @param verbose Print informative messages. Defaults to TRUE.
#' @param threshold Specify a threshold for low ACF. NULL estimates the threshold automatically.
#' @describeIn ar_acf Autocorrelation checker for `eeg_ICA` objects
#' @importFrom stats sd cor
#' @importFrom graphics abline par
#' @export
ar_acf.eeg_ICA <- function(data,
                           ms = 20,
                           plot = TRUE,
                           verbose = TRUE,
                           threshold = NULL,
                           ...) {

  time_lag <- round(data$srate * (ms / 1000))
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
#' @param data `eeg_ICA` object
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
#' @examples
#' demo_sobi <- run_ICA(demo_epochs, pca = 10)
#' ar_chanfoc(demo_sobi)
#' @return A character vector of component names that break the threshold.
#' @export
ar_chanfoc <- function(data,
                       plot = TRUE,
                       threshold = NULL,
                       verbose = TRUE,
                       measure = "max",
                       ...) {

  if (!inherits(data, "eeg_ICA")) {
    stop("Only eeg_ICA objects supported.")
  }
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
#' Detect components that load heavily on a small number of trials. Looks for
#' components that have one particular trial that has a particularly high
#' z-score.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data `eeg_ICA` object
#' @param plot Produce plot showing max z-scores and threshold for all ICA
#'   components.
#' @param threshold Specify a threshold (z-score) for high focality. NULL
#'   estimates the threshold automatically.
#' @param verbose Print informative messages.
#' @references Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide
#'   to the selection of independent components of the electroencephalogram for
#'   artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
#'   10.1016/j.jneumeth.2015.02.025
#' @examples
#' demo_sobi <- run_ICA(demo_epochs, pca = 10)
#' ar_trialfoc(demo_sobi)
#' @return A character vector of component names that break the threshold.
#' @export
ar_trialfoc <- function(data,
                        plot = TRUE,
                        threshold = NULL,
                        verbose = TRUE) {

  if (!inherits(data, "eeg_ICA")) {
    stop("Only eeg_ICA objects supported.")
  }

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
    graphics::plot(matrixStats::colMaxs(zmat))
    graphics::abline(h = threshold)
  }

  trial_foc <- channel_names(data)[matrixStats::colMaxs(zmat) > threshold]

  if (verbose) {
    message("Components with high trial focality: ", paste0(trial_foc, sep = " "))
  }

  trial_foc
}
