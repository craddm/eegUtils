#'Plot 1-D timecourse data.
#'
#'Typically event-related potentials/fields, but could also be timecourses from
#'frequency analyses for single frequencies. Averages over all submitted
#'electrodes. Output is a ggplot2 object.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @examples
#' plot_timecourse(demo_epochs, "A29")
#' plot_timecourse(demo_epochs, "A29", add_CI = TRUE)
#' @param data EEG dataset. Should have multiple timepoints.
#' @param ... Other arguments passed to methods.
#' @importFrom dplyr summarise group_by ungroup
#' @import ggplot2
#' @return Returns a ggplot2 plot object
#' @export
plot_timecourse <- function(data,
                            ...) {
  UseMethod("plot_timecourse", data)
}

#' @export
plot_timecourse.default <- function(data,
                                    ...) {
  stop("plot_timecourse() doesn't handle objects of class ",
       class(data))
}

#'@param electrode Electrode(s) to plot.
#'@param add_CI Add confidence intervals to the graph. Defaults to 95 percent
#'  between-subject CIs.
#'@param time_lim Character vector. Numbers in whatever time unit is used
#'  specifying beginning and end of time-range to plot. e.g. c(-.1, .3)
#'@param facet Deprecated. Please use standard ggplot2 facetting functions.
#'@param baseline Character vector. Times to use as a baseline. Takes the mean
#'  over the specified period and subtracts. e.g. c(-.1,0)
#'@param colour Variable to colour lines by. If no variable is passed, only one
#'  line is drawn.
#'@param color Alias for colour.
#'@describeIn plot_timecourse Plot a data.frame timecourse
#'@export
plot_timecourse.data.frame <- function(data,
                               electrode = NULL,
                               time_lim = NULL,
                               facet,
                               add_CI = FALSE,
                               baseline = NULL,
                               colour = NULL,
                               color = NULL,
                               ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  if (!is.null(electrode)) {
    data <- select_elecs(data,
                         electrode)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        baseline)
  }

  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim)
  }

  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- as.name(color)
    }
  } else {
    colour <- as.name(colour)
  }

  tc_plot <- create_tc(data,
                       add_CI = FALSE,
                       colour = colour)
  tc_plot
}

#' @describeIn plot_timecourse plot \code{eeg_evoked} timecourses
#' @export
plot_timecourse.eeg_evoked <- function(data,
                               electrode = NULL,
                               time_lim = NULL,
                               facet,
                               add_CI = FALSE,
                               baseline = NULL,
                               colour = NULL,
                               color = NULL,
                               ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  if (add_CI) {
    warning("Cannot add_CI for eeg_evoked objects.")
    add_CI <- FALSE
  }

  data <- parse_for_tc(data,
                       time_lim,
                       electrode,
                       baseline,
                       add_CI)

  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- as.name(color)
    }
  } else {
    colour <- as.name(colour)
  }

  tc_plot <- create_tc(data,
                       add_CI = add_CI,
                       colour = colour)

  tc_plot
}

#' @describeIn plot_timecourse Plot individual components from \code{eeg_ICA} components
#' @param component name or number of ICA component to plot
#' @export
plot_timecourse.eeg_ICA <- function(data,
                            component = NULL,
                            time_lim = NULL,
                            facet,
                            add_CI = FALSE,
                            baseline = NULL,
                            colour = NULL,
                            color = NULL,
                            ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  # Select specifed times
  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim = time_lim)
  }

  ## Select specified electrodes -----
  if (!is.null(component)) {
    data <- select_elecs(data,
                         component = component)
  }

  ## check for US spelling of colour...
  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- as.name(color)
    }
  } else {
    colour <- as.name(colour)
  }

  ## Do baseline correction
  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline)
  }

  if (!add_CI) {
    data <- eeg_average(data)
  }

  data <- as.data.frame(data,
                        long = TRUE)

  tc_plot <- create_tc(data,
                       add_CI = add_CI,
                       colour = colour)

  tc_plot
  }

#' @describeIn plot_timecourse Plot timecourses from \code{eeg_epochs} objects.
#' @export
plot_timecourse.eeg_epochs <- function(data,
                               electrode = NULL,
                               time_lim = NULL,
                               facet,
                               add_CI = FALSE,
                               baseline = NULL,
                               colour = NULL,
                               color = NULL, ...) {

  if (!missing(facet)) {
    warning("The facet parameter is deprecated. Please use facet_wrap/facet_grid")
    facet <- NULL
  }

  data <- parse_for_tc(data,
                       time_lim = time_lim,
                       electrode = electrode,
                       baseline = baseline,
                       add_CI = add_CI)
  ## check for US spelling of colour...
  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- as.name(color)
    }
  } else {
    colour <- as.name(colour)
  }

  tc_plot <- create_tc(data,
                       add_CI = add_CI,
                       colour = colour)

  tc_plot
}

#' @describeIn plot_timecourse Plot timecourses from \code{eeg_epochs} objects.
#'
plot_timecourse.eeg_stats <- function(data,
                                      time_lim = NULL,
                                      electrode = NULL,
                                      ...) {
  
  data <- parse_for_tc(data,
                       time_lim = time_lim,
                       electrode = electrode,
                       baseline = NULL,
                       add_CI = FALSE)
  ## check for US spelling of colour...
  # if (is.null(colour)) {
  #   if (!is.null(color)) {
  #     colour <- as.name(color)
  #   }
  # } else {
  #   colour <- as.name(colour)
  # }

  tc_plot <- create_tc(data,
                       add_CI = FALSE,
                       colour = NULL,
                       quantity = statistic)

  tc_plot
}

#' Parse data for timecourses
#'
#' Internal command for parsing various data structures into a suitable format
#' for \code{tc_plot}
#'
#' @param data data to be parsed
#' @param time_lim time limits to be returned.
#' @param electrode electrodes to be selected
#' @param baseline baseline times to be average and subtracted
#' @param add_CI Logical for whether CIS are required
#' @keywords internal
parse_for_tc <- function(data,
                         time_lim,
                         electrode,
                         baseline,
                         add_CI) {

  if (is.eeg_ICA(data) && is.null(electrode)) {
    stop("Component number must be supplied for ICA.")
  }

  ## Select specified electrodes -----
  #if (is.eeg_stats(data)) {
   # data <- select(data, electrode)
  if (!is.null(electrode)) {
    data <- select(data, electrode)
    #data <- select_elecs(data,
     #                    electrode)
  }

  ## Do baseline correction
  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline)
  }

  # Select specifed times
  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim = time_lim)
  }

  if (!is.eeg_stats(data) && !is.eeg_evoked(data) && !add_CI) {
    data <- eeg_average(data)
  }

  data <- as.data.frame(data,
                        long = TRUE)
}

#' Internal function for creation of timecourse plots
#'
#' @param data A data frame to be plotted
#' @param add_CI whether to add confidence intervals
#' @param colour whether to use colour
#' @keywords internal
create_tc <- function(data,
                      add_CI,
                      colour,
                      quantity = amplitude) {

  if (is.null(colour)) {
    tc_plot <- ggplot2::ggplot(data,
                               aes(x = time,
                                   y = {{quantity}}))
  } else {
    colour <- ggplot2::enquo(colour)
    tc_plot <- ggplot2::ggplot(data,
                               aes(x = time,
                                  y = amplitude,
                                  colour = !!colour))
  }

  if (add_CI) {
    if (is.null(colour)) {
      tc_plot <- tc_plot +
        stat_summary(fun.data = mean_cl_normal,
                     geom = "ribbon",
                     linetype = "dashed",
                     fill = NA,
                     colour = "black",
                     size = 1,
                     alpha = 0.5)
    } else {
      tc_plot <- tc_plot +
        stat_summary(fun.data = mean_cl_normal,
                     geom = "ribbon",
                     linetype = "dashed",
                     aes(colour = !!colour),
                     fill = NA,
                     size = 1,
                     alpha = 0.5)
    }
  }

  tc_plot <- tc_plot +
    stat_summary(fun.y = "mean",
                 geom = "line",
                 size = 1.2)
  tc_plot +
    labs(x = "Time (s)",
         y = expression(paste("Amplitude (", mu, "V)")),
         colour = "",
         fill = "") +
    geom_vline(xintercept = 0,
               linetype = "solid", size = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size = .5)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
}

