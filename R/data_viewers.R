#'Browse continuous data.
#'
#'A Shiny gadget for browsing continuous data interactively. Data can be viewed
#'as a butterfly plot (all electrodes overlaid) or as individual traces
#'(electrodes "stacked").
#'
#'@author Matt Craddock \email{matt@mattcraddock.com}
#'
#'@param df Data to be plotted. Must have columns headed time, amplitude, and
#'  electrode.
#'@param sig_length Length of signal to be plotted initially (seconds if
#'  continuous, epochs if epoched).
#'@param n_elecs Number of electrodes to be plotted on a single screen.
#'@param downsample Only works on \code{eeg_data} objects. Reduces size of data
#'  by only plotting every 4th point, speeding up plotting considerably.
#'
#'@import tidyr
#'@import dplyr
#'@import ggforce
#'@import ggplot2
#'@import shiny
#'@import miniUI
#'

browse_data <- function(data, sig_length = 5, n_elecs = 20, downsample = FALSE) {

  if (is.eeg_data(data)) {
    #data <- cbind(data$signals, data$timings)
    #data <- gather(data, electrode, amplitude, -sample, -time)
    data <- as.data.frame(data, long = TRUE)
    if (downsample) {
      data <- data[seq(1, nrow(data), 4), ]
    }
  }

  ui <- miniPage(
    gadgetTitleBar("Data browser"),
    miniTabstripPanel(
      miniTabPanel(title = "Butterfly", icon = icon("line-chart"),
                   miniContentPanel(
                     fillCol(
                       flex = c(4, NA, 1),
                       plotOutput("butterfly", height = '100%'),
                       sliderInput("time_range",
                                   label = "Display start time",
                                   step = 1,
                                   min = 0,
                                   max = max(unique(data$time)),
                                   value = min(unique(data$time)),
                                   width = '100%'),
                       fillRow(
                         numericInput("sig_time", "Display length", sig_length, min = 1, max = 60),
                         numericInput("uV_scale", "Scale (microvolts)", max(data$signals), min = 10),
                         checkboxInput("dc_offset", "Remove DC offset", value = TRUE)
                         )
                       )
                     )
                   ),
      miniTabPanel(title = "Individual" ,
                   miniContentPanel(
                     wellPanel(
                       plotOutput("time_plot"),
                       style = "overflow-y:scroll; max-height: 800px"
                     ),
                     fillPage(
                       fillCol(
                         flex = c(NA, 2),

                         sliderInput("time_range_ind",
                                     label = "Display start time",
                                     step = 1,
                                     min = 0,
                                     max = max(unique(data$time)),
                                     value = min(unique(data$time)),
                                     width = '100%'),
                         fillRow(
                           numericInput("sig_time_ind", "Display length", sig_length, min = 1, max = 60),
                           numericInput("elecs_per_page_ind", "Electrodes per page", n_elecs, min = 1, max = 30),
                           checkboxInput("dc_offset", "Remove DC offset", value = TRUE)
                           )
                         )
                       )
                     )
                   )
              )
    )


  server <- function(input, output, session) {

    output$butterfly <- renderPlot({
      tmp_data <- dplyr::filter(data,
                                time >= input$time_range,
                                time <= (input$time_range + input$sig_time))
      if (input$dc_offset) {
        tmp_data <- rm_baseline(tmp_data)
      }

      zz <- plot_butterfly(tmp_data, legend = FALSE, continuous = TRUE)
      zz
    })

    output$time_plot <- renderPlot({
      tmp_data <- dplyr::filter(data,
                                time >= input$time_range_ind,
                                time <= (input$time_range_ind + input$sig_time))
      if (input$dc_offset) {
        tmp_data <- rm_baseline(tmp_data)
      }

      init_plot <- ggplot2::ggplot(tmp_data, aes(x = time, y = amplitude), environment = environment()) +
        geom_line() +
        facet_grid(electrode~., scales = "free_y", switch = "y") +
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
          strip.text.y = element_text(angle = 180),
          panel.spacing = unit(0, "lines"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()
        ) +
        scale_x_continuous(expand = c(0,0))

      init_plot
    }, height = 2000, inline = FALSE)

    observeEvent(input$done, {
      stopApp()
    })
    session$onSessionEnded(stopApp)
  }
  runGadget(ui, server, viewer = paneViewer(minHeight = 600))
}
