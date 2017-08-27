#' Browse continuous data.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param df Data to be plotted.
#' @param sig_length Length of signal to be plotted initially (seconds if continuous, epochs if epoched).

#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @import ggforce
#' @import shiny
#' @import miniUI
#'

browse_data <- function(data, sig_length = 5) {

  ui <- miniPage(
    gadgetTitleBar("Data browser"),
    miniContentPanel(
      fillCol(
        flex = c(2, NA, NA),
        plotOutput("time_plot", height = "100%"),
        sliderInput("time_range",
                    label = "Display start time",
                    step = 1,
                    min = 0,
                    max = max(unique(data$time)),
                    value = min(unique(data$time)),
                    width = '100%'),
        numericInput("sig_time", "Display length", sig_length, min = 1, max = 60)
      )
    )
  )

  server <- function(input, output, session) {
    output$time_plot <- renderPlot({
      #time_ <- c(0, sig_length)
      tmp_data <- dplyr::filter(data, time >= input$time_range, time <= (input$time_range + input$sig_time))
      init_plot <- ggplot(tmp_data, aes(x = time, y = amplitude)) +
        geom_line() +
        facet_grid_paginate(electrode~.,
                            scales = "free_y",
                            switch = "y",
                            ncol = 1,
                            nrow = 20,
                            page = 1,
                            byrow = TRUE) +
        theme_bw() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_text(angle = 180)
        ) +
        scale_x_continuous(expand = c(0,0))
      init_plot
    }, height = 'auto')

    observeEvent(input$done, {
      stopApp()
    })


  }
  runGadget(ui, server, viewer = paneViewer(minHeight = 600))
}
