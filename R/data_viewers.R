#'Browse EEG data.
#'
#'A Shiny gadget for browsing continuous data interactively. Data can be viewed
#'as a butterfly plot (all electrodes overlaid) or as individual traces
#'(electrodes "stacked"). Currently, the scale cannot be manually set and is
#'determined by the range of the viewable data.
#'
#'@author Matt Craddock \email{matt@mattcraddock.com}
#'
#'@param data \code{eeg_data} object to be plotted.
#'@param sig_length Length of signal to be plotted initially (seconds if
#'  continuous, epochs if epoched).
#'@param n_elecs Number of electrodes to be plotted on a single screen. (not yet implemented)
#'@param downsample Only works on \code{eeg_data} objects. Reduces size of data
#'  by only plotting every 4th point, speeding up plotting considerably. Nb. may
#'  cause aliasing as does not filter. Use with caution.
#'
#'@importFrom dplyr filter
#'@import ggplot2
#'@import shiny
#'@import miniUI
#'@export
#'

browse_data <- function(data, sig_length = 5, n_elecs = NULL, downsample = FALSE) {

  if (!is.eeg_data(data)) {
    stop("Function only implemented for eeg_data objects.")
  } else {

    continuous <- data$continuous
    srate <- data$srate
    data <- as.data.frame(data, long = TRUE)

    if (downsample) {
      if (continuous) {
        data <- group_by(data, electrode)
        data <- mutate(data, amplitude = iir_filt(amplitude, high_freq = 0.8 * (srate / 4 / 2), srate = srate))
        data <- ungroup(data)
      } else {
        data <- group_by(data, electrode, epoch)
        data <- mutate(data, amplitude = iir_filt(amplitude, high_freq = 0.8 * (srate / 4 / 2), srate = srate))
        data <- ungroup(data)
      }
      data <- data[seq(1, nrow(data), 4), ]
    }
  }

  if (continuous) {
    ui <- miniPage(
      gadgetTitleBar("Continous data browser"),
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
                           numericInput("sig_time", "Display length", value = sig_length, min = 1, max = 60),
                           #numericInput("uV_scale", "Scale (microvolts)", value = 50, min = 1),
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
                             #numericInput("elecs_per_page_ind", "Electrodes per page", n_elecs, min = 1, max = 30),
                             checkboxInput("dc_offset_ind", "Remove DC offset", value = TRUE)
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
                                  time <= (input$time_range_ind + input$sig_time_ind))
        if (input$dc_offset_ind) {
          tmp_data <- rm_baseline(tmp_data)
        }

        init_plot <- ggplot2::ggplot(tmp_data, aes(x = time, y = amplitude)) +
          geom_line() +
          facet_grid(electrode~., scales = "free_y", switch = "y") +
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
      }, height = 2000)

      observeEvent(input$done, {
        stopApp()
      })
      session$onSessionEnded(stopApp)
    }
  } else {

    ui <- miniPage(
      gadgetTitleBar("Epoched data browser"),
      miniTabstripPanel(
        miniTabPanel(title = "Butterfly", icon = icon("line-chart"),
                     miniContentPanel(
                       fillCol(
                         flex = c(4, NA, 1),
                         plotOutput("butterfly", height = '100%'),
                         sliderInput("time_range",
                                       label = "Display starting epoch",
                                       step = 1,
                                       min = 1,
                                       max = max(unique(data$epoch)),
                                       value = min(unique(data$epoch)),
                                       width = '100%'),
                         fillRow(
                           numericInput("sig_time", "No. epochs", sig_length, min = 1, max = 60),
                           #numericInput("uV_scale", "Scale (microvolts)", max(data$amplitude), min = 10, value = 100),
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
                                       label = "Display starting epoch",
                                       step = 1,
                                       min = 1,
                                       max = max(unique(data$epoch)),
                                       value = min(unique(data$epoch)),
                                       width = '100%'),
                           fillRow(
                             numericInput("sig_time_ind", "Display length", sig_length, min = 1, max = 60),
                             #numericInput("elecs_per_page_ind", "Electrodes per page", n_elecs, min = 1, max = 30),
                             checkboxInput("dc_offset_ind", "Remove DC offset", value = TRUE)
                           )
                         )
                       )
                     )
        )
      )
    )

    server <- function(input, output, session) {

      tmp_dat <- reactive({
        dplyr::filter(data,
                      epoch >= input$time_range,
                      epoch <= (input$time_range + (input$sig_time - 1)))
      })

      tmp_data <- debounce(tmp_dat, 1000)

      output$butterfly <- renderPlot({
        #tmp_data <- dplyr::filter(data,
         #                         epoch >= input$time_range,
          #                        epoch <= (input$time_range + (input$sig_time - 1)))
        if (input$dc_offset) {
          real_data <- rm_baseline(tmp_data())
        } else {
          real_data <- tmp_data()
        }

        butter_out <- plot_butterfly(real_data, legend = FALSE, browse_mode = TRUE) +
          facet_wrap("epoch", nrow = 1) +
          theme(
            panel.spacing = unit(0, "lines")
            ) +
          geom_vline(xintercept = max(unique(real_data$time))) +
          geom_vline(xintercept = 0, linetype = "longdash")
        butter_out
      })

      tmp_dat_ind <- reactive({
        dplyr::filter(data,
                      epoch >= input$time_range_ind,
                      epoch <= (input$time_range_ind + (input$sig_time_ind - 1)))
      })

      tmp_data_ind <- debounce(tmp_dat_ind, 1000)

      output$time_plot <- renderPlot({

        if (input$dc_offset_ind) {
          real_data <- rm_baseline(tmp_data_ind())
        } else {
          real_data <- tmp_data_ind()
        }

        init_plot <- ggplot2::ggplot(real_data, aes(x = time, y = amplitude)) +
          geom_line() +
          facet_grid(electrode ~ epoch, scales = "free_y", switch = "y") +
          theme_minimal() +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            strip.text.y = element_text(angle = 180),
            panel.spacing = unit(0, "lines"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
          ) +
          scale_x_continuous(expand = c(0,0)) +
          geom_vline(xintercept = max(unique(real_data$time))) +
          geom_vline(xintercept = 0, linetype = "longdash")

        init_plot
      }, height = 2500)

      observeEvent(input$done, {
        stopApp()
      })
      session$onSessionEnded(stopApp)
    }
  }
  runGadget(ui, server, viewer = paneViewer(minHeight = 600))
}



#' Interactive 3d plot of electrode locations
#'
#' @import plotly
#'

plot_electrodes <- function (data) {
  data <- electrodeLocs(data)
  xyz <- pol_to_sph(data$theta, data$radius)
  plotly::plot_ly(x, y, z, unique(data$electrode))
}
