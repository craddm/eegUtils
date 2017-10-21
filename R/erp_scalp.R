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
                      size = .8,
                      show_guide = TRUE) {

  if (is.null(color)) {
    data <- group_by(data, electrode, time)
    data <- summarise(data, amplitude = mean(amplitude))
  } else {
    data <- group_by_(data, electrode, time, as.name(color))
    data <- summarise(data, amplitude = mean(amplitude))
  }

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
          axis.ticks = element_line(size=.3),
          plot.margin = unit(c(8,8,8,8), "pt"))
  if (show_guide) {
    p <- p + annotation_custom(grob = ggplotGrob(guide),
                               xmin = min(data$x)-.07, xmax = min(data$x)+.09,
                               ymin = min(data$y)-.09, ymax = min(data$y)+.09)
  }
  for (i in 1:nrow(data)) {
    p <- p + annotation_custom(grob = ggplotGrob(data$plot[[i]]),
                               xmin = data$x[i]-.05, xmax = data$x[i]+.05,
                               ymin = data$y[i]-.07, ymax = data$y[i]+.07)
  }

  return(p)
}

#' Interactive scalp maps
#'
#' Launches a Shiny Gadget for an interactive version of erp_scalp, allowing clicking of individual electrodes to plot them in a separate panel. In that panel, they can be averaged over or plotted as individual electrodes.
#' Currently buggy with datasets with existing electrode co-ordinates.
#'
#' @param data An EEG dataset.
#' @param colour Variable to colour lines by. If no variable is passed, only
#' one line is drawn for each electrode.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#'
#' @import shiny
#' @import miniUI
#' @export

interactive_scalp <- function(data, colour = NULL) {

  ui <- miniPage(
    gadgetTitleBar("Scalp ERPs"),
    miniTabstripPanel(
      miniTabPanel("Whole scalp", icon = icon("circle"),
                   miniContentPanel(
                     fillCol(
                       flex = c(7,1),
                       plotOutput("Scalp", height = "100%",
                                  click = "click_plot"),
                       verbatimTextOutput("click_info")
                     )
                   ),
                   miniButtonBlock(
                     actionButton("reset", "Reset selection")
                   )
      ),
      miniTabPanel("Selected electrodes", icon = icon("line-chart"),
                   miniContentPanel(
                     plotOutput("Selected", height = "100%")
                   ),
                   miniButtonBlock(
                     actionButton("avg", "Mean"),
                     actionButton("single", "Plot individual electrodes")
                   )
      )
    )
  )

  server <- function(input, output, session) {

    data <- electrode_locations(data, drop = TRUE)

    button_reacts <- reactiveValues(sel_elecs = list(), avg = TRUE)

    output$Scalp <- renderPlot({
      if (is.null(colour)) {
        erp_scalp(data)
      } else {
        erp_scalp(data, color = as.name(colour))
      }
    })

    observeEvent(input$click_plot, {

      tmp <- nearPoints(data, input$click_plot, "x", "y", threshold = 45, maxpoints = 1)

      if (nrow(tmp) > 0) {
        if (tmp$electrode %in% button_reacts$sel_elecs) {
          button_reacts$sel_elecs <- button_reacts$sel_elecs[-which(button_reacts$sel_elecs == tmp$electrode)]
        } else {
          button_reacts$sel_elecs <- c(button_reacts$sel_elecs, tmp$electrode)
        }
      }

      output$click_info <- renderPrint({
        cat("Selected:", unlist(button_reacts$sel_elecs))
      })

      # plot selected electrodes when on appropriate tab
      output$Selected <- renderPlot({
        if (button_reacts$avg){
          if (is.null(colour)) {
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ])
          } else {
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ], colour = as.name(colour))
          }
        } else {
          if (is.null(colour)) {
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ], colour = "electrode")
          } else{
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ], colour = as.name(colour)) +
              facet_wrap(~electrode)
          }
        }
      })
    })

    #Close app gracefully if done button clicked

    observeEvent(input$done, {
      stopApp()
    })

    # Check if reset button clicked on Whole scalp page
    observeEvent(input$reset, {
      button_reacts$sel_elecs <- NULL
    })

    observeEvent(input$avg, {
      button_reacts$avg <- TRUE
    })

    observeEvent(input$single, {
      button_reacts$avg <- FALSE
    })

  }

  runGadget(ui, server)
}
