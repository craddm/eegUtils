#'Browse EEG data.
#'
#'A Shiny gadget for browsing EEG data and ICA decompositions interactively.
#'With EEG data (epoched or continuous), data can be viewed as a butterfly plot
#'(all electrodes overlaid) or as individual traces (electrodes "stacked").
#'Currently, the scale cannot be manually set and is determined by the range of
#'the viewable data. With
#'
#'@author Matt Craddock \email{matt@@mattcraddock.com}
#'@import ggplot2
#'@import shiny
#'@import miniUI
#'@param data \code{eeg_data}, \code{eeg_epochs}, or \code{eeg_ICA} object to be
#'  plotted.
#'@param ... Other parameters passed to browsing functions.
#'@export

browse_data <- function(data, ...) {
  UseMethod("browse_data", data)
}

#' @export
#' @describeIn browse_data View \code{eeg_ICA} component properties
browse_data.eeg_ICA <- function(data,
                                ...) {

  ui <- miniPage(
    gadgetTitleBar("ICA dashboard"),
    miniContentPanel(
      fillRow(
        fillCol(
          plotOutput("topo_ica",
                     height = "100%"),
          plotOutput("comp_psd",
                     height = "100%")),
        fillCol(plotOutput("comp_img",
                           height = "100%"),
                plotOutput("comp_tc",
                           height = "100%"),
                shiny::selectInput("icomp",
                                   "Component",
                                   names(data$signals)))
        )
      )
    )

  server <- function(input,
                     output,
                     session) {

    output$topo_ica <- renderCachedPlot({
      comp_no <- which(names(data$signals) == input$icomp)
      topoplot(data,
               component = comp_no,
               verbose = FALSE)
      },
      cacheKeyExpr = {input$icomp})

    output$comp_img <- renderCachedPlot({
      erp_image(data,
                component = input$icomp)
    },
    cacheKeyExpr = {input$icomp})

    output$comp_tc <- renderCachedPlot({
      plot_timecourse(data,
                      component = input$icomp)
    },
    cacheKeyExpr = {input$icomp})

    output$comp_psd <- renderCachedPlot({
      tmp_psd <-
        compute_psd(select(data, input$icomp),
                    n_fft = data$srate,
                    verbose = FALSE)

      tmp_psd <- dplyr::rename(tmp_psd,
                               power = 2)
      ggplot(tmp_psd,
             aes(x = frequency,
                 y = 10 * log10((power)))) +
        stat_summary(geom = "ribbon",
                     fun.data = mean_cl_normal,
                     alpha = 0.5) +
        stat_summary(geom = "line",
                     fun.y = mean) +
        coord_cartesian(xlim = c(2, 50)) +
        theme_classic() +
        labs(x = "Frequency (Hz)", y = "Power (dB)")
        },
      cacheKeyExpr = {input$icomp})

    observeEvent(input$done, {
      returnValue <- ""
      stopApp(returnValue)
    })
    session$onSessionEnded(stopApp)
  }
  runGadget(ui,
            server,
            viewer = paneViewer(minHeight = 600))
}

#' @export
browse_data.eeg_stats <- function(data, ...) {
  warning("Not currently implemented for eeg_stats objects.")
}

#'@param sig_length Length of signal to be plotted initially (seconds if
#'  continuous, epochs if epoched).
#'@param n_elecs Number of electrodes to be plotted on a single screen. (not yet
#'  implemented)
#'@param downsample Only works on \code{eeg_data} or \code{eeg_epochs} objects.
#'  Reduces size of data by only plotting every 4th point, speeding up plotting
#'  considerably. Defaults to TRUE for eeg_data, FALSE for eeg_epochs
#'
#'@export
#'@describeIn browse_data Browse continuous EEG data.

browse_data.eeg_data <- function(data,
                                 sig_length = 5,
                                 n_elecs = NULL,
                                 downsample = TRUE,
                                 ...) {

  srate <- data$srate

  if (downsample) {
    data <- eeg_downsample(data,
                           q = 4)
  }

  ui <- miniPage(
      gadgetTitleBar("Continous data browser"),
      miniTabstripPanel(
        miniTabPanel(title = "Butterfly",
                     icon = icon("line-chart"),
                     miniContentPanel(
                       fillCol(
                         flex = c(4, NA, 1),
                         plotOutput("butterfly",
                                    height = "100%"),
                         sliderInput("time_range",
                                     label = "Display start time",
                                     step = 1,
                                     min = 0,
                                     max = max(unique(data$timings$time)),
                                     value = min(unique(data$timings$time)),
                                     width = "100%"),
                         fillRow(
                           numericInput("sig_time",
                                        "Display length",
                                        value = sig_length,
                                        min = 1,
                                        max = 60),
                           #numericInput("uV_scale", "Scale (microvolts)", value
                           #= 50, min = 1),
                           checkboxInput("dc_offset",
                                         "Remove DC offset",
                                         value = TRUE)
                         )
                       )
                     )
        ),
        miniTabPanel(title = "Individual",
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
                                       max = max(unique(data$timings$time)),
                                       value = min(unique(data$timings$time)),
                                       width = "100%"),
                           fillRow(
                             numericInput("sig_time_ind",
                                          "Display length",
                                          sig_length,
                                          min = 1, max = 60),
                             #numericInput("elecs_per_page_ind", "Electrodes per
                             #page", n_elecs, min = 1, max = 30),
                             checkboxInput("dc_offset_ind",
                                           "Remove DC offset",
                                           value = TRUE)
                           )
                         )
                       )
                     )
        )
      )
    )

    server <- function(input,
                       output,
                       session) {

      output$butterfly <- renderPlot({
        # select only the time range that we want to display
        tmp_data <- select_times(data,
                                 time_lim = c(input$time_range,
                                              input$time_range + input$sig_time))

        if (input$dc_offset) {
          tmp_data <- rm_baseline(tmp_data,
                                  verbose = FALSE)
        }

        zz <- plot_butterfly(tmp_data,
                             legend = FALSE,
                             continuous = TRUE)
        zz
      })

      output$time_plot <- renderPlot({
        tmp_data <- select_times(data,
                                 time_lim = c(input$time_range_ind,
                                              input$time_range_ind + input$sig_time_ind))

        if (input$dc_offset_ind) {
          tmp_data <- rm_baseline(tmp_data,
                                  verbose = FALSE)
        }

        tmp_data <- as.data.frame(tmp_data,
                                  long = TRUE)

        init_plot <- ggplot2::ggplot(tmp_data,
                                     aes(x = time,
                                         y = amplitude)) +
          geom_line() +
          facet_grid(electrode ~ .,
                     scales = "free_y",
                     switch = "y") +
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
          scale_x_continuous(expand = c(0, 0))

        init_plot
      }, height = 2000)

      observeEvent(input$done, {
        stopApp()
      })
      session$onSessionEnded(stopApp)
    }
  runGadget(ui,
            server,
            viewer = paneViewer(minHeight = 600))
}

#'@export
#'@describeIn browse_data Browse epoched EEG data.

browse_data.eeg_epochs <- function(data,
                                   sig_length = 5,
                                   n_elecs = NULL,
                                   downsample = FALSE,
                                   ...) {

  if (downsample) {
    data <- eeg_downsample(data,
                           q = 4)
  }

  ui <- miniPage(
    gadgetTitleBar("Epoched data browser"),
    miniTabstripPanel(
      miniTabPanel(title = "Butterfly",
                   icon = icon("line-chart"),
                   miniContentPanel(
                     fillCol(
                       flex = c(4, NA, 1),
                       plotOutput("butterfly", height = "100%"),
                       sliderInput("time_range",
                                   label = "Display starting epoch",
                                   step = 1,
                                   min = 1,
                                   max = max(unique(data$timings$epoch)),
                                   value = min(unique(data$timings$epoch)),
                                   width = "100%"),
                       fillRow(
                         numericInput("sig_time",
                                      "No. epochs",
                                      sig_length,
                                      min = 1, max = 60),
                         #numericInput("uV_scale", "Scale (microvolts)",
                         #max(data$amplitude), min = 10, value = 100),
                         checkboxInput("dc_offset",
                                       "Remove DC offset",
                                       value = FALSE)
                         )
                       )
                     )
                   ),
      miniTabPanel(title = "Individual",
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
                                     max = max(unique(data$timings$epoch)),
                                     value = min(unique(data$timings$epoch)),
                                     width = "100%"),
                         fillRow(
                           numericInput("sig_time_ind",
                                        "Display length",
                                        sig_length,
                                        min = 1, max = 60),
                           #numericInput("elecs_per_page_ind", "Electrodes per
                           #page", n_elecs, min = 1, max = 30),
                           checkboxInput("dc_offset_ind",
                                         "Remove DC offset",
                                         value = FALSE)
                           )
                         )
                       )
                     )
                   )
      )
    )

    server <- function(input,
                       output,
                       session) {

      tmp_dat <- reactive({
        select_epochs(data,
                      epoch_no = seq(input$time_range,
                                     input$time_range + input$sig_time - 1))
      })

      tmp_data <- debounce(tmp_dat,
                           1000)

      output$butterfly <- renderPlot({

        if (input$dc_offset) {
          tmp_data <- rm_baseline(tmp_data(), verbose = FALSE)
        } else {
          tmp_data <- tmp_data()
        }

        tmp_data <- as.data.frame(tmp_data,
                                  long = TRUE)

        butter_out <- plot_butterfly(tmp_data,
                                     legend = FALSE,
                                     browse_mode = TRUE) +
          facet_wrap("epoch",
                     nrow = 1) +
          theme(
            panel.spacing = unit(0, "lines")
            ) +
          geom_vline(xintercept = max(unique(tmp_data$time))) +
          geom_vline(xintercept = 0,
                     linetype = "longdash")
        butter_out
      })

      tmp_dat_ind <- reactive({
        select_epochs(data,
                      epoch_no = seq(input$time_range_ind,
                                     input$time_range_ind + input$sig_time_ind - 1))
      })

      tmp_data_ind <- debounce(tmp_dat_ind, 1000)

      output$time_plot <- renderPlot({

        if (input$dc_offset_ind) {
          tmp_data_ind <- rm_baseline(tmp_data_ind(),
                                      verbose = FALSE)
        } else {
          tmp_data_ind <- tmp_data_ind()
        }

        tmp_data_ind <- as.data.frame(tmp_data_ind,
                                      long = TRUE)

        init_plot <- ggplot2::ggplot(tmp_data_ind,
                                     aes(x = time,
                                         y = amplitude)) +
          geom_line() +
          facet_grid(electrode ~ epoch,
                     scales = "free_y",
                     switch = "y") +
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
          scale_x_continuous(expand = c(0, 0)) +
          geom_vline(xintercept = max(unique(tmp_data_ind$time))) +
          geom_vline(xintercept = 0,
                     linetype = "longdash")

        init_plot
      },
      height = 2500)

      observeEvent(input$done, {
        stopApp()
      })
      session$onSessionEnded(stopApp)
    }
  runGadget(ui,
            server,
            viewer = paneViewer(minHeight = 600))
}
