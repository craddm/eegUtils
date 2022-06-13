#' Plot ERP difference waves
#'
#' Calculates the difference between the event-related potentials from two
#' conditions and plots it.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @param data `eegUtils` object. Should have multiple timepoints.
#' @param ... Other arguments passed to methods.
#' @examples
#' plot_difference(demo_spatial, conditions = "epoch_labels", electrode = "P8")
#' @return Returns a `ggplot2` plot object
#' @export
plot_difference <- function(data,
                            ...) {
  UseMethod("plot_difference", data)
}

#' @export
plot_difference.default <- function(data,
                                    ...) {
  stop("plot_difference() doesn't handle objects of class ",
       class(data))
}

#'@param electrode Electrode(s) to plot.
#'@param conditions Defaults to "epoch_labels".
#'@param time_lim Character vector. Numbers in whatever time unit is used
#'  specifying beginning and end of time-range to plot. e.g. c(-.1, .3)
#'@param baseline Character vector. Times to use as a baseline. Takes the mean
#'  over the specified period and subtracts. e.g. c(-.1,0)
#'@param colour Variable to colour lines by. If no variable is passed, only one
#'  line is drawn.
#'@param color Alias for colour.
#'@param mapping A ggplot2 `aes()` mapping.
#'@describeIn plot_difference Plot an ERP difference wave from an `eeg_epochs` object
#' @export
plot_difference.eeg_epochs <-
  function(data,
           electrode = NULL,
           time_lim = NULL,
           baseline = NULL,
           colour = NULL,
           color = NULL,
           mapping = NULL,
           conditions = "epoch_labels",
           ...) {

    data <-
      parse_for_tc(
        data,
        time_lim = time_lim,
        electrode = electrode,
        baseline = baseline,
        add_CI = FALSE,
        mapping = mapping,
        facets = conditions,
        colour = NULL
      )

    cond_levels <- unique(data[[conditions]])

    data <-
      tidyr::pivot_wider(
        data,
        id_cols = c(participant_id,
                    time,
                    electrode),
        names_from = conditions,
        values_from = amplitude,
        values_fn = mean
        )

    if (length(cond_levels) == 2) {
      data$difference <- data[[cond_levels[1]]] - data[[cond_levels[2]]]
      diff_label <- paste(cond_levels[1], "-", cond_levels[2])
    } else {
      stop("Can only currently plot differences for two levels")
    }

    data

    tc_plot <-
      create_tc(data,
                add_CI = FALSE,
                colour = colour,
                mapping = mapping,
                quantity = difference) +
      labs(subtitle = diff_label)
    tc_plot
  }
