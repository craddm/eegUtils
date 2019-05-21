#' Create a butterfly plot from timecourse data
#'
#' Typically event-related potentials/fields, but could also be timecourses from
#' frequency analyses for single frequencies. Output is a ggplot2 object. CIs
#' not possible.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data EEG dataset. Should have multiple timepoints.
#' @param ... Other parameters passed to plot_butterfly
#' @examples
#' plot_butterfly(demo_epochs)
#' plot_butterfly(demo_epochs,
#' time_lim = c(-.1, .4),
#' legend = FALSE)
#' @export

plot_butterfly <- function(data, ...) {
  UseMethod("plot_butterfly", data)
}

#' @param time_lim Character vector. Numbers in whatever time unit is used
#'   specifying beginning and end of time-range to plot. e.g. c(-.1,.3)
#' @param baseline  Character vector. Times to use as a baseline. Takes the mean
#'   over the specified period and subtracts. e.g. c(-.1, 0)
#' @param colourmap Attempt to plot using a different colourmap (from
#'   RColorBrewer). (Not yet implemented)
#' @param legend Plot legend or not.
#' @param facet Deprecated. Please use standard ggplot2 facetting functions.
#' @param continuous Is the data continuous or not (I.e. epoched)
#' @param browse_mode Custom theme for use with browse_data.
#' @return ggplot2 object showing ERPs for all electrodes overlaid on a single
#'   plot.
#' @import ggplot2
#' @importFrom dplyr group_by ungroup summarise
#' @importFrom tidyr gather
#' @describeIn plot_butterfly Default `plot_butterfly` method for data.frames, \code{eeg_data}
#' @export

plot_butterfly.default <- function(data,
                                   time_lim = NULL,
                                   baseline = NULL,
                                   colourmap = NULL,
                                   legend = TRUE,
                                   continuous = FALSE,
                                   browse_mode = FALSE,
                                   facet,
                                   ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  if (browse_mode == FALSE) {
    data <- dplyr::group_by(data,
                            time,
                            electrode)
    data <- dplyr::summarise(data,
                             amplitude = mean(amplitude))
    data <- dplyr::ungroup(data)
  }

  ## select time-range of interest -------------

  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        baseline)
  }

  #Set up basic plot -----------
  create_bf(data,
            legend = legend,
            browse_mode = browse_mode,
            continuous = FALSE)
}

#' @describeIn plot_butterfly Plot butterfly for \code{eeg_evoked} objects
#' @export
plot_butterfly.eeg_evoked <- function(data,
                                      time_lim = NULL,
                                      baseline = NULL,
                                      colourmap = NULL,
                                      legend = TRUE,
                                      continuous = FALSE,
                                      browse_mode = FALSE,
                                      facet,
                                      ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  data <- parse_for_bf(data,
                       time_lim,
                       baseline)

  create_bf(data,
            legend = legend,
            browse_mode = browse_mode,
            continuous = FALSE)

}

#' @describeIn plot_butterfly Butterfly plot for EEG statistics

plot_butterfly.eeg_stats <- function(data,
                                     time_lim = NULL,
                                     baseline = NULL,
                                     colourmap = NULL,
                                     legend = TRUE,
                                     continuous = FALSE,
                                     browse_mode = FALSE,
                                     ...) {

  data <- data.frame(data$statistic,
                     time = data$timings)
  data <- tidyr::gather(data,
                        key = "electrode",
                        value = "amplitude",
                        -time)

  plot_butterfly(data,
                 time_lim,
                 baseline,
                 colourmap,
                 legend,
                 continuous,
                 browse_mode)

}

#' @describeIn plot_butterfly Create butterfly plot for \code{eeg_data} objects
#' @export
plot_butterfly.eeg_data <- function(data,
                                    time_lim = NULL,
                                    baseline = NULL,
                                    legend = TRUE,
                                    facet,
                                    browse_mode = FALSE,
                                    ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  data <- parse_for_bf(data,
                       time_lim,
                       baseline)
  create_bf(data,
            legend = legend,
            browse_mode = browse_mode,
            continuous = TRUE)
}

#' @describeIn plot_butterfly Create butterfly plot for \code{eeg_epochs} objects
#' @export
plot_butterfly.eeg_epochs <- function(data,
                                      time_lim = NULL,
                                      baseline = NULL,
                                      legend = TRUE,
                                      facet,
                                      browse_mode = FALSE,
                                      ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  data <- eeg_average(data)
  data <- parse_for_bf(data,
                       time_lim,
                       baseline)
  create_bf(data,
            legend = legend,
            browse_mode = browse_mode,
            continuous = FALSE)
}

#' Parse data for butterfly plots
#'
#' Internal command for parsing various data structures into a suitable format
#' for \code{plot_butterfly}
#'
#' @param data data to be parsed
#' @param time_lim time limits to be returned.
#' @param baseline baseline times to be average and subtracted
#' @keywords internal
parse_for_bf <- function(data,
                         time_lim = NULL,
                         baseline = NULL) {

  # Select specifed times
  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim = time_lim)
  }

  ## Do baseline correction
  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline)
  }
  data <- as.data.frame(data,
                        long = TRUE,
                        coords = FALSE)
  data
}

#' @import ggplot2
#' @keywords internal
create_bf <- function(data,
                      legend,
                      browse_mode,
                      continuous) {

  #Set up basic plot -----------
  butterfly_plot <-
    ggplot2::ggplot(data,
                    aes(x = time,
                        y = amplitude))

  if (length(unique(data$epoch)) > 1) {
    chan_lines <- function() {
      stat_summary(fun.y = mean,
                   geom = "line",
                   aes(colour = electrode),
                   alpha = 0.5)
    }
  } else {
    chan_lines <- function() {
      geom_line(aes(colour = electrode),
                alpha = 0.5)
    }
  }


  if (browse_mode) {
    butterfly_plot <-
      butterfly_plot +
      geom_line(colour = "black",
                aes(group = electrode),
                alpha = 0.2) +
      labs(x = "Time (s)",
           y = expression(paste("Amplitude (", mu, "V)")),
           colour = "") +
      geom_hline(yintercept = 0,
                 size = 0.5,
                 linetype = "dashed",
                 alpha = 0.5) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.ticks = element_line(size = .5))
  } else {
    butterfly_plot <-
      butterfly_plot +
      # geom_line(alpha = 0.5,
      #           aes(colour = electrode)) +
      chan_lines() +
      labs(x = "Time (s)",
           y = expression(paste("Amplitude (", mu, "V)")),
           colour = "") +
      geom_hline(yintercept = 0, size = 0.5) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.ticks = element_line(size = .5))

    if (!continuous) {
      butterfly_plot <-
        butterfly_plot +
        geom_vline(xintercept = 0, size = 0.5)
    }
  }

  if (legend) {
    butterfly_plot +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))
  } else {
    butterfly_plot +
      theme(legend.position = "none")
  }
}
