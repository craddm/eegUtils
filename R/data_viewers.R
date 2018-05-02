#'Browse EEG data.
#'
#'A Shiny gadget for browsing EEG data interactively. Data can be viewed
#'as a butterfly plot (all electrodes overlaid) or as individual traces
#'(electrodes "stacked"). Currently, the scale cannot be manually set and is
#'determined by the range of the viewable data.
#'
#'@author Matt Craddock \email{matt@@mattcraddock.com}
#'
#'@param data \code{eeg_data} object to be plotted.
#'@param ... Other parameters passed to browsing functions.
#'@export

browse_data <- function(data, ...) {
  UseMethod("browse_data", data)
}

#' @export
browse_data.eeg_ICA <- function(data, ...) {
  warning("Not currently implemented for eeg_ICA objects.")
}

#' @export
browse_data.eeg_stats <- function(data, ...) {
  warning("Not currently implemented for eeg_stats objects.")
}

#'@param sig_length Length of signal to be plotted initially (seconds if
#'  continuous, epochs if epoched).
#'@param n_elecs Number of electrodes to be plotted on a single screen. (not yet
#'  implemented)
#'@param downsample Only works on \code{eeg_data} objects. Reduces size of data
#'  by only plotting every 4th point, speeding up plotting considerably.
#'
#'@import dplyr
#'@import ggplot2
#'@import shiny
#'@import miniUI
#'@export
#'@describeIn browse_data Browse continuous EEG data.

browse_data.eeg_data <- function(data, sig_length = 5, n_elecs = NULL,
                        downsample = FALSE, ...) {

  continuous <- data$continuous
  srate <- data$srate
  #data <- as.data.frame(data, long = TRUE)
  if (downsample) {
    data <- eeg_downsample(data, q = 4)
    #data <- iir_filt(data, high_freq = 0.8 * (data$srate / 4 / 2))
    #data <- as.data.frame(data, long = TRUE)
    #data <- data[seq(1, nrow(data), 4), ]
    } else {
    #data <- as.data.frame(data, long = TRUE)
    }

  ui <- miniPage(
      gadgetTitleBar("Continous data browser"),
      miniTabstripPanel(
        miniTabPanel(title = "Butterfly", icon = icon("line-chart"),
                     miniContentPanel(
                       fillCol(
                         flex = c(4, NA, 1),
                         plotOutput("butterfly", height = "100%"),
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
                                        min = 1, max = 60),
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

    server <- function(input, output, session) {

      output$butterfly <- renderPlot({
        # select only the time range that we want to display
        tmp_data <- select_times(data,
                                 time_lim = c(input$time_range,
                                              input$time_range + input$sig_time))

        if (input$dc_offset) {
          tmp_data <- rm_baseline(tmp_data)
        }

        zz <- plot_butterfly(tmp_data, legend = FALSE, continuous = TRUE)
        zz
      })

      output$time_plot <- renderPlot({
        tmp_data <- select_times(data,
                                 time_lim = c(input$time_range_ind,
                                              input$time_range_ind + input$sig_time_ind))

        if (input$dc_offset_ind) {
          tmp_data <- rm_baseline(tmp_data)
        }

        tmp_data <- as.data.frame(tmp_data, long = TRUE)

        init_plot <- ggplot2::ggplot(tmp_data, aes(x = time, y = amplitude)) +
          geom_line() +
          facet_grid(electrode ~ ., scales = "free_y", switch = "y") +
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
  runGadget(ui, server, viewer = paneViewer(minHeight = 600))
}

#'@export
#'@describeIn browse_data Browse epoched EEG data.

browse_data.eeg_epochs <- function(data, sig_length = 5,
                                   n_elecs = NULL, downsample = FALSE, ...) {

  if (downsample) {
    data <- eeg_downsample(data, q = 4)
    #data <- iir_filt(data, high_freq = 0.8 * (data$srate / 4 / 2))
    #data <- as.data.frame(data, long = TRUE)
    #data <- data[seq(1, nrow(data), 4), ]
  } else {
    #data <- as.data.frame(data, long = TRUE)
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

    server <- function(input, output, session) {

      tmp_dat <- reactive({
        select_epochs(data,
                      epoch_no = seq(input$time_range,
                                     input$time_range + input$sig_time - 1))
      })

      tmp_data <- debounce(tmp_dat, 1000)

      output$butterfly <- renderPlot({

        if (input$dc_offset) {
          tmp_data <- rm_baseline(tmp_data())
        } else {
          tmp_data <- tmp_data()
        }

        tmp_data <- as.data.frame(tmp_data, long = TRUE)

        butter_out <- plot_butterfly(tmp_data, legend = FALSE,
                                     browse_mode = TRUE) +
          facet_wrap("epoch", nrow = 1) +
          theme(
            panel.spacing = unit(0, "lines")
            ) +
          geom_vline(xintercept = max(unique(tmp_data$time))) +
          geom_vline(xintercept = 0, linetype = "longdash")
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
          tmp_data_ind <- rm_baseline(tmp_data_ind())
        } else {
          tmp_data_ind <- tmp_data_ind()
        }

        tmp_data_ind <- as.data.frame(tmp_data_ind, long = TRUE)

        init_plot <- ggplot2::ggplot(tmp_data_ind,
                                     aes(x = time, y = amplitude)) +
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
          scale_x_continuous(expand = c(0, 0)) +
          geom_vline(xintercept = max(unique(tmp_data_ind$time))) +
          geom_vline(xintercept = 0, linetype = "longdash")

        init_plot
      },
      height = 2500)

      observeEvent(input$done, {
        stopApp()
      })
      session$onSessionEnded(stopApp)
    }
  runGadget(ui, server, viewer = paneViewer(minHeight = 600))
}



#' Plot electrode locations
#'
#' Produces either a 2D plot of the electrode locations or an interactive plot
#' of electrode locations in 3D space.
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @importFrom plotly plot_ly
#' @param data Data with associated electrode locations to be plotted.
#' @param interact Choose 2D cartesian layout, or, if set to TRUE, an
#'   interactive 3D plot of electrode locations. Defaults to FALSE.
#' @param ... other parameters
#' @export

plot_electrodes <- function(data, interact = FALSE, ...) {
  UseMethod("plot_electrodes", data)
}

#' @importFrom plotly plot_ly
#' @import ggplot2
#' @describeIn plot_electrodes generic plot electrodes function

plot_electrodes.default <- function(data, interact = FALSE, ...) {
  if ("electrode" %in% names(data)) {
    data <- data.frame(electrode = unique(data$electrode))
    data <- electrode_locations(data)

    if (interact) {
      xyz <- pol_to_sph(data$theta, data$radius)
      xyz$electrode <- unique(data$electrode)
      plotly::plot_ly(xyz, x = ~x, y = ~y, z = ~z, text = ~electrode,
                      type = "scatter3d", mode = "text+markers")
    } else {
      ggplot2::ggplot(data, aes(x = x, y = y, label = electrode)) +
        geom_text() +
        theme_minimal() +
        coord_equal()
    }
  } else {
    stop("No electrodes found.")
  }
}

#' @importFrom plotly plot_ly
#' @describeIn plot_electrodes Plot electrodes associated with an \code{eeg_data} object.
plot_electrodes.eeg_data <- function(data, interact = FALSE, ...) {

  if (is.null(data$chan_info)) {
    warning("Adding standard locations...")
    data <- electrode_locations(data)
  }

  if (interact) {
    xyz <- pol_to_sph(data$chan_info$theta, data$chan_info$radius)
    xyz$electrode <- data$chan_info$electrode
    plotly::plot_ly(xyz, x = ~x, y = ~y, z = ~z, text = ~electrode,
                    type = "scatter3d", mode = "text+markers")
  } else {
    ggplot2::ggplot(data$chan_info, aes(x = x, y = y, label = electrode)) +
      geom_text() +
      theme_minimal() +
      coord_equal()
  }
}
