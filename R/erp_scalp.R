#' Plot ERPs on the scalp
#'
#' Creates an ERP figure for each electrode and layouts them on the scalp.
#'
#' @param data An EEG dataset.
#' @param electrode Column name containing electrode names in data.
#' Defaults to "electrode".
#' @param amplitude Column name containing amplitudes in data.
#' Defaults to "amplitude".
#' @param time Column name containing time in data. Defaults to "time".
#' @param color Variable to color lines by. If no variable is passed, only
#' one line is drawn for each electrode.
#' @param size Size of line(s).
#' @param show_guide Should a guide be shown.
#'
#' @details The function uses default electrode names and locations contained
#' in the package.
#'
#' @author Matti Vuorre, \email{mv2521@columbia.edu}
#' @importFrom purrr map
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @return Returns a ggplot2 plot object.
#' @export

erp_scalp <- function(data,
                      electrode = "electrode",
                      amplitude = "amplitude",
                      time = "time",
                      color = NULL,
                      size = .65,
                      show_guide = TRUE) {

  data <- as.data.frame(data)

  # Data maxima for plot limits
  maxAmp <- max(data[, amplitude])
  minAmp <- min(data[, amplitude])
  maxTime <- max(data[, time])
  minTime <- min(data[, time])

  plotfun <- function(x) {
    plot <- ggplot(data = x, aes_(as.name(time), as.name(amplitude)))
    plot <- plot + facet_wrap("electrodefacet")
    plot <- plot + coord_cartesian(ylim = c(minAmp, maxAmp),
                                   xlim = c(minTime, maxTime))
    plot <- plot + geom_vline(xintercept = 0, size=.3)
    plot <- plot + geom_hline(yintercept = 0, size=.3)
    plot <- plot + theme_void()
    plot <- plot + theme(strip.text = element_text(size=8))
    if (!is.null(color)) {
      plot <- plot + scale_color_brewer(palette = "Set1")
      plot <- plot + geom_line(aes_(color = as.name(color)),
                               show.legend = FALSE, size = size)
    } else {
      plot <- plot + geom_line(size = size)
    }
    return(plot)
  }

  data$electrodefacet <- data[, electrode]
  data <- nest(data, -electrode)
  data <- mutate(data, plot = map(data, plotfun))
  data <- select(data, -data)

  # Get default electrode locations from pkg internal data
  data <- electrode_locations(data, drop = T)

  minx <- min(data$x)
  maxx <- max(data$x)
  miny <- min(data$y)
  maxy <- max(data$y)


  # Are these lines necessary? Commented out to help with interactive plotting -
  # clicks in Shiny will use native co-ordinate space, easier to then match with electrode names
  #----

  #data$x <- data$x + abs(minx)
  #data$y <- data$y + abs(miny)

  p <- ggplot(data, aes(x, y)) +
    geom_blank() +
    theme_void() +
    theme(plot.margin = unit(c(8,8,8,8), "pt"))

  guide <- ggplot(data, aes(x=time, y=amplitude)) +
    coord_cartesian(ylim = c(minAmp, maxAmp),
                    xlim = c(minTime, maxTime)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    geom_vline(xintercept = 0, size = .4) +
    geom_hline(yintercept = 0, size = .4) +
    labs(y = expression(paste("Amplitude (", mu,"V)")),
         x = "Time (ms)") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.3))
  if (show_guide) {
    p <- p + annotation_custom(grob = ggplotGrob(guide),
                               xmin = min(data$x)-.07, xmax = min(data$x)+.09,
                               ymin = min(data$y)-.07, ymax = min(data$y)+.07)
  }
  for (i in 1:nrow(data)) {
    p <- p + annotation_custom(grob = ggplotGrob(data$plot[[i]]),
                               xmin = data$x[i]-.055, xmax = data$x[i]+.055,
                               ymin = data$y[i]-.055, ymax = data$y[i]+.055)
  }

  return(p)
}
