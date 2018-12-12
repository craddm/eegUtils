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
#' @importFrom rlang parse_quo
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
#'@param baseline Character vector. Times to use as a baseline. Takes the mean
#'  over the specified period and subtracts. e.g. c(-.1,0)
#'@param colour Variable to colour lines by. If no variable is passed, only one
#'  line is drawn.
#'@param color Alias for colour.
#'@param group (not yet implemented)
#'@param facet Create multiple plots for a specified grouping variable.
#'@describeIn plot_timecourse Plot a data.frame timecourse
#'@export
plot_timecourse.data.frame <- function(data,
                               electrode = NULL,
                               time_lim = NULL,
                               group = NULL,
                               facet = NULL,
                               add_CI = FALSE,
                               baseline = NULL,
                               colour = NULL,
                               color = NULL,
                               ...) {

  if (!is.null(electrode)) {
    data <- select_elecs(data, electrode)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data, baseline)
  }

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
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

}

#' @describeIn plot_timecourse plot \code{eeg_evoked} timecourses
#' @export
plot_timecourse.eeg_evoked <- function(data,
                               electrode = NULL,
                               time_lim = NULL,
                               group = NULL,
                               facet = NULL,
                               add_CI = FALSE,
                               baseline = NULL,
                               colour = NULL,
                               color = NULL,
                               ...) {

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
                            group = NULL,
                            facet = NULL,
                            add_CI = FALSE,
                            baseline = NULL,
                            colour = NULL,
                            color = NULL,
                            ...) {

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
                        long = T)

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
                               group = NULL,
                               facet = NULL,
                               add_CI = FALSE,
                               baseline = NULL,
                               colour = NULL,
                               color = NULL, ...) {

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


#' @describeIn plot_tc plot_tc for eeg_stats objects.
#' @noRd
#'
plot_timecourse.eeg_stats <- function(data, time_lim, ...) {
  warning("Not yet implemented.")
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

  if (is.eeg_ICA(data) & is.null(electrode)) {
    stop("Component number must be supplied for ICA.")
  }

  ## Select specified electrodes -----
  if (!is.null(electrode)) {
    data <- select_elecs(data,
                         electrode)
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

  if (!is.eeg_evoked(data) & !add_CI) {
    data <- eeg_average(data)
  }

  data <- as.data.frame(data, long = TRUE)
}

#' Internal function for creation of timecourse plots
#'
#' @param data A data frame to be plotted
#' @param add_CI whether to add confidence intervals
#' @param colour whether to use colour
#' @keywords internal
create_tc <- function(data,
                      add_CI,
                      colour) {

  if (is.null(colour)) {
    tc_plot <- ggplot2::ggplot(data,
                               aes(x = time,
                                   y = amplitude))
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
    labs(x = "Time (s)", y = expression(paste("Amplitude (", mu, "V)")),
         colour = "", fill = "") +
    geom_vline(xintercept = 0, linetype = "solid", size = 0.5) +
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

#' Create a butterfly plot from timecourse data
#'
#' Typically event-related potentials/fields, but could also be timecourses from
#' frequency analyses for single frequencies. Output is a ggplot2 object. CIs
#' not possible.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data EEG dataset. Should have multiple timepoints.
#' @param ... Other parameters passed to plot_butterfly
#' @export

plot_butterfly <- function(data, ...) {
  UseMethod("plot_butterfly", data)
}

#' @param time_lim Character vector. Numbers in whatever time unit is used
#'   specifying beginning and end of time-range to plot. e.g. c(-.1,.3)
#' @param group Group lines by a specificed grouping variable.
#' @param baseline  Character vector. Times to use as a baseline. Takes the mean
#'   over the specified period and subtracts. e.g. c(-.1, 0)
#' @param facet Create multiple plots for a specified grouping variable.
#' @param colourmap Attempt to plot using a different colourmap (from
#'   RColorBrewer). (Not yet implemented)
#' @param legend Plot legend or not.
#' @param continuous Is the data continuous or not (I.e. epoched)
#' @param browse_mode Custom theme for use with browse_data.
#' @return ggplot2 object showing ERPs for all electrodes overlaid on a single
#'   plot.
#' @import ggplot2
#' @importFrom dplyr group_by ungroup summarise
#' @import tidyr
#' @describeIn plot_butterfly Default `plot_butterfly` method for data.frames, \code{eeg_data}
#' @export

plot_butterfly.default <- function(data,
                           time_lim = NULL,
                           group = NULL,
                           facet = NULL,
                           baseline = NULL,
                           colourmap = NULL,
                           legend = TRUE,
                           continuous = FALSE,
                           browse_mode = FALSE,
                           ...) {

  if (browse_mode == FALSE && is.null(facet)) {
    data <- dplyr::group_by(data, time, electrode)
    data <- dplyr::summarise(data, amplitude = mean(amplitude))
    data <- dplyr::ungroup(data)
  }

  ## select time-range of interest -------------

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data, baseline)
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
                                      group = NULL,
                                      facet = NULL,
                                      baseline = NULL,
                                      colourmap = NULL,
                                      legend = TRUE,
                                      continuous = FALSE,
                                      browse_mode = FALSE,
                                      ...) {

  if (identical(class(data$signals), "list")) {
    time_vec <- data$timings$time
    data <- Reduce("+", data$signals) / length(data$signals)
    data$time <- time_vec
    data <- tidyr::gather(data,
                          electrode,
                          amplitude,
                          -time,
                          factor_key = T)
  } else {
    data <- as.data.frame(data, long = TRUE)
  }

  plot_butterfly(data,
                 time_lim,
                 group,
                 facet,
                 baseline,
                 colourmap,
                 legend,
                 continuous,
                 browse_mode)
  }

#' @describeIn plot_butterfly Butterfly plot for EEG statistics

plot_butterfly.eeg_stats <- function(data,
                                     time_lim = NULL,
                                     group = NULL,
                                     facet = NULL,
                                     baseline = NULL,
                                     colourmap = NULL,
                                     legend = TRUE,
                                     continuous = FALSE,
                                     browse_mode = FALSE,
                                     ...) {

  data <- data.frame(data$statistic,
                     time = data$timings)
  data <- tidyr::gather(data, key = "electrode", value = "amplitude", -time)

  plot_butterfly(data, time_lim,
                 group,
                 facet,
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
                             browse_mode = FALSE,
                             ...) {

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
                               browse_mode = FALSE,
                               ...) {

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
                        long = TRUE)
  data
}

#' @import ggplot2
#' @keywords internal
create_bf <- function(data,
                      legend,
                      browse_mode,
                      facet = FALSE,
                      continuous) {

  #Set up basic plot -----------
  butterfly_plot <- ggplot2::ggplot(data,
                                    aes(x = time,
                                        y = amplitude))

  if (browse_mode) {
    butterfly_plot <- butterfly_plot +
      geom_line(aes(group = electrode),
                colour = "black",
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
    butterfly_plot <- butterfly_plot +
      geom_line(aes(group = electrode,
                    colour = electrode),
                alpha = 0.5) +
      labs(x = "Time (s)",
           y = expression(paste("Amplitude (", mu, "V)")),
           colour = "") +
      geom_hline(yintercept = 0, size = 0.5) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.ticks = element_line(size = .5))

    if (!continuous) {
      butterfly_plot <- butterfly_plot +
        geom_vline(xintercept = 0, size = 0.5)
    }

    if (!is.null(facet)) {
      butterfly_plot +
        facet_wrap(facet)
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
