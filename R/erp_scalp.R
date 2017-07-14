#' Plot ERPs on the scalp
#'
#' Creates an ERP figure for each electrode and layouts them on the scalp.
#'
#' @param d An EEG dataset.
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
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @return Returns a ggplot2 plot object.
#' @export

erp_scalp <- function(d,
                      electrode = "electrode",
                      amplitude = "amplitude",
                      time = "time",
                      color = NA,
                      size = .65,
                      show_guide = TRUE) {

  d <- as.data.frame(d)

  # Data maxima for plot limits
  maxAmp <- max(d[, amplitude])
  minAmp <- min(d[, amplitude])
  maxTime <- max(d[, time])
  minTime <- min(d[, time])

  plotfun <- function(x) {
    plot <- ggplot(data = x, aes_(as.name(time), as.name(amplitude)))
    plot <- plot + facet_wrap("electrodefacet")
    plot <- plot + coord_cartesian(ylim = c(minAmp, maxAmp),
                                   xlim = c(minTime, maxTime))
    plot <- plot + geom_vline(xintercept = 0, size=.3)
    plot <- plot + geom_hline(yintercept = 0, size=.3)
    plot <- plot + theme_void()
    plot <- plot + theme(strip.text = element_text(size=8))
    if (!is.na(color)) {
      plot <- plot + scale_color_brewer(palette = "Set1")
      plot <- plot + geom_line(aes_(color = as.name(color)),
                               show.legend = FALSE, size = size)
    } else {
      plot <- plot + geom_line(size = size)
    }
    return(plot)
  }

  d$electrodefacet <- d[, electrode]
  d <- nest(d, -electrode)
  d <- mutate(d, plot = map(data, plotfun))
  d <- select(d, -data)

  # Get default electrode locations from pkg internal data
  d <- join_locs(d, .drop = T)

  minx <- min(d$x)
  maxx <- max(d$x)
  miny <- min(d$y)
  maxy <- max(d$y)

  d$x <- d$x + abs(minx)
  d$y <- d$y + abs(miny)

  p <- ggplot(d, aes(x, y)) +
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
                               xmin = min(d$x)-.07, xmax = min(d$x)+.09,
                               ymin = min(d$y)-.07, ymax = min(d$y)+.07)
  }
  for (i in 1:nrow(d)) {
    p <- p + annotation_custom(grob = ggplotGrob(d$plot[[i]]),
                               xmin = d$x[i]-.055, xmax = d$x[i]+.055,
                               ymin = d$y[i]-.055, ymax = d$y[i]+.055)
  }

  return(p)
}
