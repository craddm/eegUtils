#' Plot 1d timecourse data.
#'
#' Typically event-related potentials/fields, but could also be timecourses from frequency analyses for single frequencies. Averages over all submitted electrodes. Output is a ggplot2 object.
#'
#'@author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#'
#'@param data EEG dataset. Should have multiple timepoints.
#'@param electrode Electrode(s) to plot.
#'@param add_CI Add confidence intervals to the graph. Defaults to 95 percent between-subject CIs.
#'@param time_lim Character vector. Numbers in whatever time unit is used specifying beginning and end of time-range to plot. e.g. c(-100,300)
#'@param baseline Character vector. Times to use as a baseline. Takes the mean over the specified period and subtracts. e.g. c(-100,0)
#'@param colour Variable to colour lines by. If no variable is passed, only
#' one line is drawn for each electrode.
#'@param facet Create multiple plots for a specified grouping variable.
#'
#'@import dplyr
#'@import ggplot2
#'@return Returns a ggplot2 plot object
#'@export

plot_timecourse <- function(data, time_lim = NULL, group = NULL, facet = NULL, add_CI = FALSE, baseline = NULL, colour = NULL, color = NULL, electrode = NULL) {

  # Filter out unwanted timepoints, and find nearest time values in the data --------------

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  # Average over all epochs in data (respecting "conditions"). --
  if ("epoch" %in% colnames(data)){
    data <- summarise(group_by(data, time, electrode, condition),
                      amplitude = mean(amplitude))
  }

  if (!is.null(group)) {
    if (group %in% colnames(data)){

    }
  }

  ## Set up basic plot -----------
  tc_plot <- ggplot(data, aes(x = time, y = amplitude))

  if (!(is.null(colour) & is.null(color))) {
    tc_plot <- tc_plot + scale_color_brewer(palette = "Set1")
    tc_plot <- tc_plot +
      stat_summary(fun.y = mean,
                   geom = "line",
                   aes_(colour = as.name(colour)),
                   size = 1)
  } else {
    tc_plot <- tc_plot +
      stat_summary(fun.y = mean,
                   geom = "line",
                   size = 1)
  }

  if (!is.null(facet)) {
    if (facet %in% colnames(data)){
      data <- mutate_(data, f1 = facet)
      tc_plot <- tc_plot %+% data + facet_wrap(~f1)
    } else {
      warning("Unrecognised column name.")
    }
  }

  ## Draw confidence intervals on plot.

  if (add_CI) {
    if (!(is.null(colour) & is.null(color))){
    tc_plot <- tc_plot +
      stat_summary(fun.data = mean_cl_normal,
                   geom = "ribbon",
                   linetype = "dashed",
                   aes_(colour = as.name(colour)),
                   fill = NA,
                   size = 1,
                   alpha = 0.5)
    } else {
      tc_plot <- tc_plot +
        stat_summary(fun.data = mean_cl_normal,
                   geom = "ribbon",
                   linetype = "dashed",
                   fill = NA,
                   size = 1,
                   alpha = 0.5)
    }
  }

  tc_plot <- tc_plot +
    labs(x = "Time (ms)", y = expression(paste("Amplitude (", mu, "V)")), colour = "", fill = "") +
    geom_vline(xintercept = 0, linetype = "solid", size = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0,0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0,0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.5))

  return(tc_plot)
}

#' Create a butterfly plot from timecourse data
#'
#' Typically event-related potentials/fields, but could also be timecourses from frequency analyses for single frequencies. Output is a ggplot2 object. CIs not possible.
#'
#' @param data EEG dataset. Should have multiple timepoints
#' @param time_lim Character vector. Numbers in whatever time unit is used specifying beginning and end of time-range to plot. e.g. c(-100,300)
#'@param baseline  Character vector. Times to use as a baseline. Takes the mean over the specified period and subtracts. e.g. c(-100,0)
#'@param facet Create multiple plots for a specified grouping variable.
#'
#' @export

plot_butterfly <- function(data, time_lim = NULL, group = NULL, facet = NULL, baseline = NULL, colourmap = NULL) {

  ## select time-range of interest -------------

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  #Set up basic plot -----------
  butterfly_plot <- ggplot(data, aes(x = time, y = amplitude))+
    geom_line(aes(group = electrode, colour = electrode))

  butterfly_plot +
    labs(x = "Time (ms)", y = expression(paste("Amplitude (", mu, "V)")), colour = "") +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_hline(yintercept = 0, size = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.5))
}
