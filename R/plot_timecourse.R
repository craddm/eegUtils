#' Plot 1d timecourse data.
#'
#' Typically event-related potentials/fields, but could also be timecourses from frequency analyses for single frequencies. Averages over all submitted electrodes. Output is a ggplot2 object.
#'
#'@author Matt Craddock, \email{matt@mattcraddock.com}
#'
#'@param data EEG dataset. Should have multiple timepoints.
#'@param electrode Electrode(s) to plot.
#'@param add_CI Add confidence intervals to the graph. Defaults to 95 percent between-subject CIs.
#'@param time_lim Character vector. Numbers in whatever time unit is used specifying beginning and end of time-range to plot. e.g. c(-100,300)
#'@param baseline Character vector. Times to use as a baseline. Takes the mean over the specified period and subtracts. e.g. c(-100,0)
#'@param colour Variable to colour lines by. If no variable is passed, only
#' one line is drawn.
#'@param group (not yet implemented)
#'@param facet Create multiple plots for a specified grouping variable.
#'
#'@import dplyr
#'@import ggplot2
#'@importFrom rlang parse_quosure
#'
#'@return Returns a ggplot2 plot object
#'
#'@export

plot_timecourse <- function(data, time_lim = NULL,
                            group = NULL, facet = NULL,
                            add_CI = FALSE, baseline = NULL,
                            colour = NULL, electrode = NULL,
                            color = NULL) {

  if (is.eeg_data(data)) {
    data <- as.data.frame(data, long = TRUE)
  }

  ## Select specified electrodes -----
  if (!is.null(electrode)) {
    data <- select_elecs(data, electrode)
  }

  ## check for US spelling of colour...
  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- color
      tmp_col <- rlang::parse_quosure(colour)
    }
  } else {
    tmp_col <- rlang::parse_quosure(colour)
  }

  ## Filter out unwanted timepoints, and find nearest time values in the data -----

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  ## Average over all epochs in data (respecting "conditions"). --
  if (is.null(colour)) {
    data <- dplyr::summarise(dplyr::group_by(data, time, electrode),
                      amplitude = mean(amplitude))
  } else {
    if ("electrode" %in% c(colour, color)) {
      data <- dplyr::summarise(dplyr::group_by(data, time, electrode),
                        amplitude = mean(amplitude))
      } else {
        data <- dplyr::summarise(dplyr::group_by(data, time, electrode, !!tmp_col),
                        amplitude = mean(amplitude))
      }
  }

  if (!is.null(group)) {
    if (group %in% colnames(data)) {

    }
  }

  ## Set up basic plot -----------
  tc_plot <- ggplot2::ggplot(data, aes(x = time, y = amplitude))

  if (!(is.null(colour) & is.null(color))) {

    tc_plot <- tc_plot + scale_color_brewer(palette = "Set1")
    tc_plot <- tc_plot +
      stat_summary(fun.y = mean,
                   geom = "line",
                   aes_(colour = as.name(colour)),
                   size = 1.2)
  } else {
    tc_plot <- tc_plot +
      stat_summary(fun.y = mean,
                   geom = "line",
                   size = 1.2)
  }

  if (!is.null(facet)) {
    if (facet %in% colnames(data)){
      data <- dplyr::mutate_(data, f1 = facet)
      tc_plot <- tc_plot %+% data + facet_wrap(~f1)
    } else {
      warning("Unrecognised column name.")
    }
  }

  ## Draw confidence intervals on plot.

  if (add_CI) {
    if (is.null(colour)) {
      tc_plot <- tc_plot +
        stat_summary(fun.data = mean_cl_normal,
                     geom = "ribbon",
                     linetype = "dashed",
                     fill = NA,
                     size = 1,
                     alpha = 0.5)
      } else {
        tc_plot <- tc_plot +
          stat_summary(fun.data = mean_cl_normal,
                       geom = "ribbon",
                       linetype = "dashed",
                       aes_(colour = as.name(colour)),
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
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = c(0,0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.5))

  return(tc_plot)
}

#' Create a butterfly plot from timecourse data
#'
#' Typically event-related potentials/fields, but could also be timecourses from frequency analyses for single frequencies. Output is a ggplot2 object. CIs not possible.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @param data EEG dataset. Should have multiple timepoints.
#' @param time_lim Character vector. Numbers in whatever time unit is used specifying beginning and end of time-range to plot. e.g. c(-100,300)
#' @param baseline  Character vector. Times to use as a baseline. Takes the mean over the specified period and subtracts. e.g. c(-100,0)
#' @param facet Create multiple plots for a specified grouping variable. (Not yet implemented)
#' @param colourmap Attempt to plot using a different colourmap (from RColorBrewer). (Not yet implemented)
#' @param legend Plot legend or not.
#' @param continuous Is the data continuous or not (I.e. epoched)
#' @return ggplot2 object showing ERPs for all electrodes overlaid on a single plot.
#' @import ggplot2
#'
#' @export

plot_butterfly <- function(data,
                           time_lim = NULL,
                           group = NULL,
                           facet = NULL,
                           baseline = NULL,
                           colourmap = NULL,
                           legend = TRUE,
                           continuous = FALSE) {

  ## select time-range of interest -------------

  if (!is.null(time_lim)) {
    data <- select_times(data, time_lim)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data, baseline)
  }

  #Set up basic plot -----------
  butterfly_plot <- ggplot(data, aes(x = time, y = amplitude))+
    geom_line(aes(group = electrode, colour = electrode), alpha = 0.5)

  butterfly_plot <- butterfly_plot +
    labs(x = "Time (ms)", y = expression(paste("Amplitude (", mu, "V)")), colour = "") +
    geom_hline(yintercept = 0, size = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.5))

  if (!continuous) {
    butterfly_plot <- butterfly_plot + geom_vline(xintercept = 0, size = 0.5)
  }

  if (legend) {
    butterfly_plot
  } else {
    butterfly_plot +
      theme(legend.position = "none")
  }
}
