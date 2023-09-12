#' Plot Global Field Power of EEG Signals
#'
#' Global Field Power (Lehmann & Skrandies, 1980) is a way to quantify the
#' amount of activity in the overall electroencephalographic signal at a given
#' timepoint. It corresponds to the spatial standard deviation.
#'
#' @param data An `eeg_epochs` object
#' @param cols Condition columns from the `epochs` metadata to calculate GFP
#'   separately for different conditions.
#' @param keep_trials Calculate GFP for each epoch separately, then average over
#'   epochs.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @export
#' @examples
#' plot_gfp(demo_spatial)
#' plot_gfp(demo_spatial, keep_trials = TRUE)
#' plot_gfp(demo_spatial, cols = "epoch_labels")
#' plot_gfp(demo_spatial, cols = "epoch_labels", keep_trials = TRUE)
plot_gfp <- function(data,
                     cols = NULL,
                     keep_trials = FALSE) {
  UseMethod("plot_gfp", data)
}

#' @export
plot_gfp.default <- function(data,
                             cols = NULL,
                             keep_trials = FALSE) {
  message("Only `eeg_epochs` objects are supported.")
}

#' @export
plot_gfp.eeg_epochs <- function(data,
                                cols = NULL,
                                keep_trials = FALSE) {

  if (keep_trials) {
    df <- as.data.frame(data)
  } else {
    if (is.null(cols)) {
      df <- eeg_average(data,
                        cols = "participant_id")
    } else {
      df <- eeg_average(data,
                        cols = c("participant_id", cols)
                        )
    }
    df <- as.data.frame(df)
  }

  chan_names <- channel_names(data)

  gfp <- cbind(
    df[, c("time", cols), drop = FALSE],
    gfp = matrixStats::rowSds(as.matrix(df[, chan_names]))
    )

  if (is.character(cols)) {
    cols <- rlang::sym(cols)
  }
  ggplot2::ggplot(gfp,
                  ggplot2::aes(x = time,
                               y = gfp,
                               colour = {{cols}})) +
    ggplot2::stat_summary(geom = "line",
                          fun.data = mean_se) +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(y = "Global field power (microvolts)",
         x = "Time (s)") +
    ggplot2::coord_cartesian(expand = FALSE)
}
