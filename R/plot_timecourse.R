#' Plot 1d timecourse data.
#'
#' Typically event-related potentials/fields, but could also be timecourses from frequency analyses for single frequencies. Averages over all submitted electrodes. Output is a ggplot2 object.
#'
#'@param df EEG dataset. Should have multiple timepoints.
#'@param add_CI Add confidence intervals to the graph. Defaults to 95 percent between-subject CIs.
#'@param time_lim Character vector. Numbers in whatever time unit is used specifying beginning and end of time-range to plot. e.g. c(-100,300)
#'@param baseline  Character vector. Times to use as a baseline. Takes the mean over the specified period and subtracts. e.g. c(-100,0)
#'@param facet Create multiple plots for a specified grouping variable.
#'
#'@import dplyr
#'@import ggplot2
#'@return Returns a ggplot2 plot object
#'@export

plot_timecourse <- function(df, time_lim = NULL, group = NULL, facet = NULL, add_CI = FALSE, baseline = NULL) {

  # Filter out unwanted timepoints, and find nearest time values in the data --------------

  if ("time" %in% colnames(df)) {
    if (length(time_lim) == 1) {
      warning("Must enter two timepoints when specifying a time range to plot; plotting whole range.")
    } else if (length(time_lim) == 2) {
      time_lim[1] <- df$time[which.min(abs(df$time - time_lim[1]))]
      time_lim[2] <- df$time[which.min(abs(df$time - time_lim[2]))]
      df <- filter(df, time >= time_lim[1] & time <= time_lim[2])
    } else {
    }
  }

  if ("epoch" %in% colnames(df)){
    df <- summarise(group_by(df, time, electrode, condition), amplitude = mean(amplitude))
  }

  if (!is.null(group)) {
    if (group %in% colnames(df)){

    }
  }


  #Set up basic plot -----------
  tc_plot <- ggplot(df, aes(x = time, y = amplitude))+
    scale_color_brewer(palette = "Set1")

  ## Draw confidence intervals on plot.



  tc_plot <- tc_plot +
    stat_summary(fun.y = mean, geom = "line", size = 1)

  if (!is.null(facet)) {
    if (facet %in% colnames(df)){
      df <- df %>%
        mutate_(f1 = facet)
      tc_plot <- tc_plot %+% df + facet_wrap(~f1)
    } else {
      warning("Unrecognised column name.")
    }
  }

  if (add_CI) {
    tc_plot <- tc_plot +
      stat_summary(fun.data = mean_cl_normal, geom = "ribbon", linetype = "dotted", size = 1, alpha = 0.5)
  }

  tc_plot +
    labs(x = "Time (ms)", y = expression(paste("Amplitude (", mu, "V)")), colour = "", fill = "") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic()


}

#' Create a butterfly plot from timecourse data
#'
#' Typically event-related potentials/fields, but could also be timecourses from frequency analyses for single frequencies. Output is a ggplot2 object. CIs not possible.
#'
#' @param df EEG dataset. Should have multiple timepoints
#' @param time_lim Character vector. Numbers in whatever time unit is used specifying beginning and end of time-range to plot. e.g. c(-100,300)
#'@param baseline  Character vector. Times to use as a baseline. Takes the mean over the specified period and subtracts. e.g. c(-100,0)
#'@param facet Create multiple plots for a specified grouping variable.
#'
#' @export

plot_butterfly <- function(df, time_lim = NULL, group = NULL, facet = NULL, baseline = NULL, colourmap = NULL) {

  ## select time-range of interest -------------

  if ("time" %in% colnames(df)) {
    if (length(time_lim) == 1) {
      warning("Must enter two timepoints when specifying a time range to plot; plotting whole range.")
    } else if (length(time_lim) == 2) {
      time_lim[1] <- df$time[which.min(abs(df$time - time_lim[1]))]
      time_lim[2] <- df$time[which.min(abs(df$time - time_lim[2]))]
      df <- filter(df, time >= time_lim[1] & time <= time_lim[2])
    } else {
    }
  }

  #Set up basic plot -----------
  butterfly_plot <- ggplot(df, aes(x = time, y = amplitude))+
    geom_line(aes(group = electrode, colour = electrode))

  butterfly_plot +
    labs(x = "Time (ms)", y = expression(paste("Amplitude (", mu, "V)")), colour = "") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic()
}
