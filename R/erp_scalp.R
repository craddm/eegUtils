#' Plot ERPs on the scalp
#'
#' Creates an ERP figure for each channel and layouts them on the scalp.
#'
#' @param d An EEG dataset.
#' @param electrode Column name containing electrode names in data.
#' @param amplitude Column name containing amplitudes in data.
#' @param time Column name containing time in data.
#' @param color Variable to color lines by.
#' @param size Size of line(s).
#' @param show_guide Should a guide be shown.
#'
#' @author Matti Vuorre, \email{mv2521@columbia.edu}
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @export

erp_scalp <- function(d,
                      electrode = "electrode",
                      amplitude = "amplitude",
                      time = "time",
                      color = "color",
                      size = .7,
                      show_guide = TRUE) {

  # Get default electrode locations from pkg internal data
  d <- as.data.frame(d)
  d <- join_locs(d, .drop = TRUE)

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
    plot <- plot + scale_color_brewer(palette = "Set1")
    if (!is.na(color)) {
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
  d <- join_locs(d, .drop = T)
  d <- mutate(d, x = x + abs(min(x)), y = y + abs(min(y)))

  p <- ggplot(d, aes(x, y)) + theme_void()

  guide <- ggplot(data, aes(x=time, y=amplitude)) +
    coord_cartesian(ylim = c(minAmp, maxAmp),
                    xlim = c(minTime, maxTime)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.5))
  if (show_guide) {
    p <- p + annotation_custom(grob = ggplotGrob(guide),
                               xmin = min(d$x)-.05, xmax = min(d$x)+.15,
                               ymin = min(d$y)-.05, ymax = min(d$y)+.1)
  }
  for (i in 1:59) {
    p <- p + annotation_custom(grob = ggplotGrob(d$plot[[i]]),
                               xmin = d$x[i]-.04, xmax = d$x[i]+.04,
                               ymin = d$y[i]-.04, ymax = d$y[i]+.04)
  }

  return(p)
}
