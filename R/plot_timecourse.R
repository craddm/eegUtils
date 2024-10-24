#'Plot one-dimensional timecourse data.
#'
#'Typically event-related potentials/fields, but could also be timecourses from
#'frequency analyses for single frequencies. Averages over all submitted
#'electrodes. For group data, `plot_timecourse` will average within-participants
#'first, using weighted averaging where possible, then across participants using
#'unweighted averaging. Output is a `ggplot2` object.
#'
#'@author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @examples
#' library(ggplot2)
#' plot_timecourse(demo_epochs, "A29")
#' plot_timecourse(demo_epochs, "A29", baseline = c(-.1, 0))
#' plot_timecourse(demo_epochs, "A29", baseline = c(-.1, 0), add_CI = TRUE)
#' plot_timecourse(demo_spatial, "Oz", baseline = c(-.1, 0), mapping = aes(colour = epoch_labels))
#' plot_timecourse(demo_spatial, "Oz", baseline = c(-.1, 0), facets = ~epoch_labels)
#'@param data EEG dataset. Should have multiple timepoints.
#'@param ... Other arguments passed to methods.
#'@import ggplot2
#'@return Returns a ggplot2 plot object
#'@export
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
#'@param mapping A `ggplot2` `aes()` mapping.
#'@param facets A right-hand-side only formula specifying which variables should
#'  be used to create facets.
#'@describeIn plot_timecourse Plot a data.frame timecourse
#'@export
plot_timecourse.data.frame <- function(data,
                                       electrode = NULL,
                                       time_lim = NULL,
                                       add_CI = FALSE,
                                       baseline = NULL,
                                       colour = NULL,
                                       color = NULL,
                                       mapping = NULL,
                                       facets = NULL,
                                       ...) {

  if (!is.null(electrode)) {
    data <- select_elecs(data,
                         electrode)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline)
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
                       colour = colour,
                       mapping = mapping,
                       facets = facets)
  tc_plot
}

#' @describeIn plot_timecourse plot `eeg_evoked` timecourses
#' @export
plot_timecourse.eeg_evoked <- function(data,
                                       electrode = NULL,
                                       time_lim = NULL,
                                       add_CI = FALSE,
                                       baseline = NULL,
                                       colour = NULL,
                                       color = NULL,
                                       mapping = NULL,
                                       facets = NULL,
                                       ...) {

  if (add_CI) {
    warning("Cannot add_CI for eeg_evoked objects.")
    add_CI <- FALSE
  }

  data <- parse_for_tc(data,
                       time_lim = time_lim,
                       electrode = electrode,
                       baseline = baseline,
                       add_CI = add_CI,
                       facets = facets,
                       mapping = mapping,
                       colour = colour)

  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- as.name(color)
    }
  } else {
    colour <- as.name(colour)
  }

  tc_plot <- create_tc(data,
                       add_CI = add_CI,
                       colour = colour,
                       mapping = mapping,
                       facets = facets)

  tc_plot
}

#' @describeIn plot_timecourse Plot individual components from `eeg_ICA`
#'   components
#' @param component name or number of ICA component to plot
#' @export
plot_timecourse.eeg_ICA <- function(data,
                                    component = NULL,
                                    time_lim = NULL,
                                    add_CI = FALSE,
                                    baseline = NULL,
                                    colour = NULL,
                                    color = NULL,
                                    mapping = NULL,
                                    facets = NULL,
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
                        long = TRUE,
                        coords = FALSE)

  tc_plot <- create_tc(data,
                       add_CI = add_CI,
                       colour = colour,
                       mapping = mapping,
                       facets = facets)

  tc_plot
}

#' @describeIn plot_timecourse Plot timecourses from `eeg_epochs` objects.
#' @export
plot_timecourse.eeg_epochs <- function(data,
                                       electrode = NULL,
                                       time_lim = NULL,
                                       add_CI = FALSE,
                                       baseline = NULL,
                                       colour = NULL,
                                       color = NULL,
                                       mapping = NULL,
                                       facets = NULL,
                                       ...) {

  data <- parse_for_tc(data,
                       time_lim = time_lim,
                       electrode = electrode,
                       baseline = baseline,
                       add_CI = add_CI,
                       facets = facets,
                       mapping = mapping,
                       colour = colour)

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
                       colour = colour,
                       mapping = mapping,
                       facets = facets)

  tc_plot
}

#' @describeIn plot_timecourse Plot timecourses from `eeg_group` objects.
#' @export
plot_timecourse.eeg_group <- function(data,
                                      electrode = NULL,
                                      time_lim = NULL,
                                      add_CI = FALSE,
                                      baseline = NULL,
                                      colour = NULL,
                                      color = NULL,
                                      mapping = NULL,
                                      facets = NULL,
                                      ...) {

  if (inherits(data,
               "eeg_tfr")) {
    return(plot_timecourse.eeg_tfr(data,
                                   electrode = electrode,
                                   time_lim = time_lim,
                                   add_CI = add_CI,
                                   baseline = baseline,
                                   colour = colour,
                                   color = color,
                                   mapping = mapping,
                                   facets = facets,
                                   ...))
  }
  data <- parse_for_tc(data,
                       time_lim = time_lim,
                       electrode = electrode,
                       baseline = baseline,
                       add_CI = add_CI,
                       facets = facets,
                       mapping = mapping,
                       colour = colour)

  if (is.null(colour)) {
    if (!is.null(color)) {
      colour <- as.name(color)
    }
  } else {
    colour <- as.name(colour)
  }

  tc_plot <- create_tc(data,
                       add_CI = add_CI,
                       colour = colour,
                       mapping = mapping,
                       facets = facets)

  tc_plot
}

#' @describeIn plot_timecourse Plot timecourses from `eeg_tfr` objects.
#' @param freq_range Choose a specific frequency range to plot. If NULL,
#'   calculates the mean over all frequencies. Note that this does not imply
#'   that there is power at an included frequency. For example, lower
#'   frequencies will have shorter timecourses than high frequencies.
#' @param type Type of baseline correction to use for `eeg_tfr` objects
#' @export
plot_timecourse.eeg_tfr <- function(data,
                                    electrode = NULL,
                                    time_lim = NULL,
                                    add_CI = FALSE,
                                    baseline = NULL,
                                    colour = NULL,
                                    color = NULL,
                                    mapping = NULL,
                                    freq_range = NULL,
                                    type = "divide",
                                    ...) {

  if (!is.null(colour) || !is.null(color)) {
    warning(
      "colour argument is kept for compatability, please use the `mapping` argument and supply a `ggplot2` `aes()` mapping. the colour parameter will be deprectated in v0.9.0 of eegUtils"
    )
  }

  if (add_CI) {
    message("Confidence intervals are not currently supported for `eeg_tfr` objects.")
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline,
                        type = type)
    yintercept <-
      switch(type,
             divide = 1,
             db = 0,
             absolute = 0,
             pc = 0,
             ratio = 1)
    ylabel <-
      switch(type,
             divide = "Power ratio",
             db = "Decibels (dB)",
             ratio = "Power ratio",
             absolute = "Power (a.u.)",
             pc = "Percent change (%)")
  } else {
    yintercept <- 0
    ylabel <- "Power (a.u.)"
  }

  if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim)
  }

  if (!is.null(electrode)) {
    data <- select_elecs(data,
                         electrode)
  }

  if (!is.null(freq_range)) {
    data <- select_freqs(data,
                         freq_range)
  }

  data_f <- as.data.frame(data,
                          long = TRUE,
                          coords = FALSE)

  tc_plot <-
    ggplot(data_f,
           aes(x = time,
               y = power)) +
    stat_summary(geom = "line",
                 fun = mean) +
    labs(x = "Time (s)",
         y = ylabel,
         colour = "",
         fill = "") +
    geom_vline(xintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    geom_hline(yintercept = yintercept,
               linetype = "solid",
               linewidth = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(linewidth = .5)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))


  if (!is.null(mapping)) {
    tc_plot <-
      tc_plot +
      mapping
  }
  tc_plot
}

#' Parse data for timecourses
#'
#' Internal command for parsing various data structures into a suitable format
#' for `tc_plot`
#'
#' @param data data to be parsed
#' @param time_lim time limits to be returned.
#' @param electrode electrodes to be selected
#' @param baseline baseline times to be average and subtracted
#' @param add_CI Logical for whether CIS are required
#' @param facets A RHS-only formula for use with `ggplot2::facet_wrap`
#' @param colour A character vector indicating which variable to use for colour.
#' @param mapping A `ggplot2` `aes()` call with axis mappings
#' @keywords internal
parse_for_tc <- function(data,
                         time_lim,
                         electrode,
                         baseline,
                         add_CI,
                         facets,
                         mapping,
                         colour) {

  if (is.eeg_ICA(data) && is.null(electrode)) {
    stop("Component number must be supplied for ICA.")
  }

  col_names <- NULL

  ## Select specified electrodes -----
  if (!is.null(electrode)) {
    data <- select(data,
                   dplyr::all_of(electrode))
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

  if (!is.null(mapping)) {
    col_names <- unname(
      vapply(mapping,
             rlang::as_label,
             character(1))
    )
  }

  if (!is.null(colour)) {
    warning(
      "colour argument is kept for compatability, please use the `mapping` argument and supply a `ggplot2` `aes()` mapping. the colour parameter will be deprectated in v0.9.0 of eegUtils"
    )
    col_names <- c(col_names, colour)
  }

  if (is.character(facets)) facets <- stats::reformulate(facets)

  col_names <- c(col_names, all.vars(facets))

  if (!is.eeg_stats(data) && !add_CI) {
    data <- eeg_average(data,
                        cols = col_names)
  }

  data <- as.data.frame(data,
                        long = TRUE,
                        coords = FALSE)
}

#' Internal function for creation of timecourse plots
#'
#' @param data A data frame to be plotted
#' @param add_CI whether to add confidence intervals
#' @param colour whether to use colour
#' @param quantity The name of the column/quantity to plot
#' @param mapping A ggplot2 `aes()` mapping.
#' @keywords internal
create_tc <- function(data,
                      add_CI,
                      colour,
                      quantity = amplitude,
                      mapping = NULL,
                      facets = NULL) {

  if (is.null(colour)) {
    tc_plot <- ggplot2::ggplot(data,
                               aes(x = time,
                                   y = {{quantity}}))
  } else {
    colour <- ggplot2::enquo(colour)
    tc_plot <- ggplot2::ggplot(data,
                               aes(x = time,
                                   y = {{quantity}},
                                   colour = {{colour}}))
  }

  if (add_CI) {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop(
        "Package \"Hmisc\" must be installed to add confidence intervals.",
        call. = FALSE
      )
    }
    if (is.null(mapping)) {
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
                       linewidth = 1,
                       alpha = 0.5)
      }
    } else {
      tc_plot <-
        tc_plot +
        stat_summary(fun.data = mean_cl_normal,
                     geom = "ribbon",
                     mapping = mapping,
                     fill = NA,
                     linewidth = 1,
                     alpha = 0.5)
    }
  }

  tc_plot <-
    tc_plot +
    stat_summary(fun = "mean",
                 geom = "line",
                 linewidth = 1.2)

  if (!is.null(mapping)) {
    tc_plot <-
      tc_plot +
      mapping
  }

  if (!is.null(facets)) {
    tc_plot <-
      tc_plot +
      facet_wrap(facets)
  }

  tc_plot +
    labs(x = "Time (s)",
         y = expression(paste("Amplitude (", mu, "V)")),
         colour = "",
         fill = "") +
    geom_vline(xintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    geom_hline(yintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(linewidth = .5)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
}
