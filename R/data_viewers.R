#' Browse continuous data.
#'
#' @author Matt Craddock \email{m.p.craddock@leeds.ac.uk}
#' @param df Data to be plotted. Must have columns headed time, amplitude, and electrode.
#' @param sig_length Length of signal to be plotted initially (seconds if continuous, epochs if epoched).
#' @param n_elecs Number of electrodes to be plotted on a single screen.
#'
#' @import tidyr
#' @import dplyr
#' @import ggforce
#' @import ggplot2
#' @import shiny
#' @import miniUI
#'

browse_data <- function(data, sig_length = 5, n_elecs = 20) {

  if (is.eeg_data(data)) {
    data <- cbind(data$signals, data$timings)
    data <- gather(data, electrode, amplitude, -sample, -time)
  }

  ui <- miniPage(
    gadgetTitleBar("Data browser"),
    miniContentPanel(
      fillCol(
        flex = c(4, NA, 1, 1),
        plotOutput("time_plot", height = "100%"),
        sliderInput("time_range",
                    label = "Display start time",
                    step = 1,
                    min = 0,
                    max = max(unique(data$time)),
                    value = min(unique(data$time)),
                    width = '100%'),
        fillRow(
          numericInput("sig_time", "Display length", sig_length, min = 1, max = 60),
          numericInput("elec_page", "Page", 1, min = 1, max = ceiling(length(unique(data$electrode))/n_elecs)),
          numericInput("elecs_per_page", "Electrodes per page", n_elecs, min = 1, max = 30)
          ),
        fillRow(
          checkboxInput("dc_offset", "Remove DC offset", value = FALSE)
          )
        )
      )
    )


  server <- function(input, output, session) {
    output$time_plot <- renderPlot({
      tmp_data <- dplyr::filter(data, time >= input$time_range, time <= (input$time_range + input$sig_time))
      if (input$dc_offset) {
        tmp_data <- rm_baseline(tmp_data)
      }
      init_plot <- ggplot(tmp_data, aes(x = time, y = amplitude), environment = environment()) +
        geom_line() +
        facet_grid(electrode~.,
                   scales = "free_y",
                   switch = "y") +
        # ggforce::facet_grid_paginate(electrode~.,
        #                       scales = "free_y",
        #                       switch = "y",
        #                       ncol = 1,
        #                       nrow = input$elecs_per_page,
        #                       page = input$elec_page,
        #                       byrow = TRUE) +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_text(angle = 180)
        ) +
        scale_x_continuous(expand = c(0,0))

      print(init_plot)
    })

    observeEvent(input$done, {
      stopApp()
    })

  }
  runGadget(ui, server, viewer = paneViewer(minHeight = 600))
}
