#' Plot ERPs on the scalp
#'
#' Creates an ERP figure for each channel and layouts them on the scalp.
#'
#' @param data An EEG dataset.
#' @param electrode The column name containing electrode names in data.
#'
#' @import dplyr
#' @export

erp_scalp <- function(data) {

  # Get default electrode locations from pkg internal data
  data <- join_locs(data, .drop = TRUE)

  # Data maxima for plot limits
  maxAmp <- max(data$amplitude)
  minAmp <- min(data$amplitude)
  maxTime <- max(data$time)
  minTime <- min(data$time)

  plotfun <- function(x) {
    ggplot(data = x, aes(time, amplitude)) +
      geom_line() +
      facet_wrap("electrodelabel") +
      coord_cartesian(ylim = c(minAmp, maxAmp),
                      xlim = c(minTime, maxTime)) +
      geom_vline(xintercept = 0, size=.3) +
      geom_hline(yintercept = 0, size=.3) +
      theme_void()
  }

  d <- data %>%
    mutate(electrodelabel = electrode) %>%
    nest(-electrode) %>%
    mutate(plot = map(data, plotfun)) %>%
    select(-data) %>%
    join_locs(.drop = T)
  d <- mutate(d, x = x + abs(min(x)), y = y + abs(min(y)))

  p <- ggplot(d, aes(x, y)) + theme_void()
  p + geom_point() + geom_label(aes(label = electrode))
  library(scales)
  guide <- ggplot(data, aes(x=time, y=amplitude)) +
    coord_cartesian(ylim = c(minAmp, maxAmp),
                    xlim = c(minTime, maxTime)) +
    scale_x_continuous(breaks = pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = pretty_breaks(n = 4)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size=.5))
  guide
  p <- p + annotation_custom(grob = ggplotGrob(guide),
                             xmin = min(d$x)-.05, xmax = min(d$x)+.15,
                             ymin = min(d$y)-.05, ymax = min(d$y)+.1)
  for (i in 1:59) {
    p <- p + annotation_custom(grob = ggplotGrob(d$plot[[i]]),
                               xmin = d$x[i]-.04, xmax = d$x[i]+.04,
                               ymin = d$y[i]-.04, ymax = d$y[i]+.04)
  }

  return(p)
}
