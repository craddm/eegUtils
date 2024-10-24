#'Browse EEG data
#'
#'A Shiny gadget for browsing EEG data and ICA decompositions interactively.
#'With EEG data (epoched or continuous), data can be viewed as a butterfly plot
#'(all electrodes overlaid) or as individual traces (electrodes "stacked"). For
#'`eeg_ICA` objects, you will instead be shown a composite of multiple
#'properties of the decomposition - a topography, an ERP image, an ERP, and a
#'power spectral density plot from 4-50 Hz.
#'
#'@author Matt Craddock \email{matt@@mattcraddock.com}
#'@import ggplot2
#'@param data `eeg_data`, `eeg_epochs`, or `eeg_ICA` object to be plotted.
#'@param ... Other parameters passed to browsing functions.
#'@export

browse_data <- function(data, ...) {
  UseMethod("browse_data", data)
}

#' @export
#' @return A character vector of component names selected for rejection.
#' @describeIn browse_data View `eeg_ICA` component properties
browse_data.eeg_ICA <- function(data,
                                ...) {

  ui <- bslib::page_fillable(
    bslib::card(
      shiny::tags$span(
        "ICA Dashboard",
        shiny::actionButton("done", "Done", class = "btn-success", style ="align-right")
        ),
      min_height = "85px"
    ),
    bslib::layout_columns(
        shiny::selectInput("icomp",
                           label = NULL,
                           names(data$signals)),
        shiny::radioButtons("reject_comps",
                            label = NULL,
                            choices = c("Keep", "Reject"),
                            inline = TRUE),
        ),
    bslib::layout_column_wrap(
      shiny::plotOutput("topo_ica", height = "400px"),
      shiny::plotOutput("comp_img"),
      shiny::plotOutput("comp_psd"),
      shiny::plotOutput("comp_tc"),
      width = 1/2,
      heights_equal = "all"
      )
  )

  server <- function(input,
                     output,
                     session) {

    comp_status <- shiny::reactiveValues()
    output$topo_ica <- shiny::bindCache(
      shiny::renderPlot({
        comp_no <- which(names(data$signals) == input$icomp)
        topoplot(data,
                 component = comp_no,
                 verbose = FALSE,
                 grid_res = 70)
        }),
        input$icomp)

    output$comp_img <- shiny::bindCache(
      shiny::renderPlot({
        erp_image(data, component = input$icomp)
        }),
        input$icomp)

    output$comp_tc <- shiny::bindCache(
      shiny::renderPlot({
        plot_timecourse(data, component = input$icomp)
        }),
      input$icomp)

    output$comp_psd <- shiny::bindCache(
      shiny::renderPlot({
        tmp_psd <-
          compute_psd(select(data, input$icomp),
                      n_fft = data$srate,
                      verbose = FALSE)
        tmp_psd <- dplyr::rename(tmp_psd,
                                 power = 2)
        tmp_psd <- dplyr::filter(tmp_psd,
                                 frequency >= 3,
                                 frequency <= 50)
        ggplot(tmp_psd,
               aes(x = frequency,
                   y = 10 * log10((power)))) +
          stat_summary(geom = "ribbon",
                       fun.data = mean_se,
                       alpha = 0.5) +
          stat_summary(geom = "line",
                       fun = mean) +
          theme_classic() +
          labs(x = "Frequency (Hz)", y = "Power (dB)")
        }),
      input$icomp)

    shiny::observeEvent(input$icomp, {
      shiny::updateRadioButtons(
        inputId = "reject_comps",
        choices = c("Keep", "Reject"),
        selected = comp_status[[input$icomp]],
        inline = TRUE)
      })

    shiny::observeEvent(input$reject_comps, {
      comp_status[[shiny::isolate(input$icomp)]] <- input$reject_comps
      comp_status
      })
    shiny::observeEvent(input$done, {
      returnValue <- shiny::reactiveValuesToList(comp_status)
      returnValue <-
        names(returnValue)[vapply(returnValue,
                                  function(x) identical(x, "Reject"),
                                  logical(1))]
      shiny::stopApp(returnValue)
    })

    session$onSessionEnded(shiny::stopApp)
  }
  shiny::runGadget(ui,
            server,
            viewer = shiny::paneViewer(minHeight = 600))
}

#' @export
browse_data.eeg_stats <- function(data, ...) {
  warning("Not currently implemented for `eeg_stats` objects.")
}

#'@param sig_length Length of signal to be plotted initially (seconds if
#'  continuous, epochs if epoched).
#'@param downsample Only works on `eeg_data` or `eeg_epochs` objects.
#'  Reduces size of data by only plotting every 4th point, speeding up plotting
#'  considerably. Defaults to TRUE for `eeg_data`, FALSE for `eeg_epochs`
#'
#'@export
#'@describeIn browse_data Browse continuous EEG data.

browse_data.eeg_data <- function(data,
                                 sig_length = 5,
                                 downsample = TRUE,
                                 ...) {

  srate <- data$srate

  if (downsample) {
    data <- eeg_downsample(data,
                           q = 4)
  }

  ui <- bslib::page_fillable(
    title = ("Continuous data browser"),
    bslib::navset_tab(
      bslib::nav_panel(
        title = "Butterfly",
        icon = shiny::icon("chart-line"),
        bslib::card(
          shiny::plotOutput("butterfly"),
          style = "resize:vertical;"
        )
      ),
      bslib::nav_panel(
        title = "Individual",
        icon = shiny::icon("chart-line"),
        bslib::card(
          shiny::wellPanel(
            shiny::plotOutput("time_plot"),
            style = "overflow-y:scroll; max-height: 800px; resize:vertical",
          ),
        )
      ),
      footer =
        bslib::card(
          shiny::sliderInput(
            "time_range",
            label = "Display time range (seconds)",
            step = 1,
            min = round(min(data$timings$time), 2),
            max = round(max(data$timings$time), 2),
            value = c(min(data$timings$time),
                      min(data$timings$time) + 4),
            width = "100%"),
          shiny::checkboxInput("dc_offset",
                               "Remove DC offset",
                               value = TRUE),
          shiny::numericInput("uV_scale",
                              "Scale (microvolts)",
                              value = 50, step = 1),
          )
      )
    )

  server <- function(input,
                     output,
                     session) {


    output$butterfly <- shiny::renderPlot({
      # select only the time range that we want to display
      tmp_data <- select_times(data,
                               time_lim = c(input$time_range[1],
                                            input$time_range[2]))

      if (input$dc_offset) {
        tmp_data <- rm_baseline(tmp_data,
                                verbose = FALSE)
      }

      butterfly_plot <- plot_butterfly(tmp_data,
                                       legend = FALSE,
                                       continuous = TRUE)
      butterfly_plot +
        coord_cartesian(ylim = c(-input$uV_scale, input$uV_scale),
                        expand = FALSE)
    })

    output$time_plot <- shiny::renderPlot({
      tmp_data <- select_times(data,
                               time_lim = c(input$time_range,
                                            input$time_range + input$sig_time))

      if (input$dc_offset) {
        tmp_data <- rm_baseline(tmp_data,
                                verbose = FALSE)
      }

      tmp_data <- as.data.frame(tmp_data,
                                long = TRUE,
                                coords = FALSE)

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
        scale_x_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(-input$uV_scale, input$uV_scale),
                        expand = FALSE)

      init_plot
    }, height = 2000)

    shiny::observeEvent(input$done, {
      stopApp()
    })
    session$onSessionEnded(stopApp)
  }
  shiny::runGadget(ui,
                   server,
                   viewer = shiny::paneViewer(minHeight = 600))
}

#'@export
#'@describeIn browse_data Browse epoched EEG data.

browse_data.eeg_epochs <- function(data,
                                   sig_length = 5,
                                   downsample = FALSE,
                                   ...) {

  if (downsample) {
    data <- eeg_downsample(data,
                           q = 4)
  }

  ui <- bslib::page_fillable(
    title = "Epoched data browser",
    bslib::navset_tab(
      bslib::nav_panel(
        title = bslib::tooltip(
          shiny::tags$span(
            "Butterfly",
            shiny::icon("chart-line")
          ),
          "Shows traces from every electrode overlapping, with transparency"),
        bslib::card(
          shiny::plotOutput("butterfly"),
          style = "resize:vertical;"
        ),
      ),
      bslib::nav_panel(
        title = bslib::tooltip(
          shiny::tags$span(
            "Individual",
            shiny::icon("chart-line"),
          ),
          "Show traces for each electrode separately"
        ),
        bslib::card(
          shiny::plotOutput("time_plot"),
          style = "overflow-y:scroll; max-height: 800px; resize:vertical;",
        )
      ),
      footer = bslib::card(
        shiny::sliderInput("time_range",
                           label = "Epoch range",
                           step = 1,
                           min = 1,
                           max = max(data$timings$epoch),
                           value = c(min(data$timings$epoch),
                                     min(data$timings$epoch) + 4),
                           width = "100%"),
        shiny::checkboxInput("dc_offset",
                             "Remove DC offset",
                             value = FALSE),
        shiny::numericInput("uV_scale",
                            "Scale (microvolts)",
                            value = 50, step = 1),
      )
    )
  )

  server <- function(input,
                     output,
                     session) {

    tmp_data <- shiny::debounce(
      shiny::reactive({
        select_epochs(data,
                      epoch_no = seq(input$time_range[1],
                                     input$time_range[2]))
      }),
      400)

    output$butterfly <- shiny::renderPlot({

      if (input$dc_offset) {
        tmp_data <- rm_baseline(tmp_data(),
                                verbose = FALSE)
      } else {
        tmp_data <- tmp_data()
      }

      tmp_data <- as.data.frame(tmp_data,
                                long = TRUE,
                                coords = FALSE)

      butter_out <- plot_butterfly(tmp_data,
                                   legend = FALSE,
                                   browse_mode = TRUE,
                                   allow_facets = TRUE) +
        facet_wrap("epoch",
                   nrow = 1) +
        theme(
          panel.spacing = unit(0, "lines")
        ) +
        geom_vline(xintercept = max(unique(tmp_data$time))) +
        geom_vline(xintercept = 0,
                   linetype = "longdash",
                   alpha = 0.2) +
        coord_cartesian(ylim = c(-input$uV_scale, input$uV_scale),
                        expand = FALSE)
      butter_out
    })

    output$time_plot <- shiny::bindCache(
      shiny::renderPlot({
        if (input$dc_offset) {
          tmp_data <- rm_baseline(tmp_data(),
                                  verbose = FALSE)
        } else {
          tmp_data <- tmp_data()
        }

        tmp_data <- as.data.frame(tmp_data,
                                  long = TRUE,
                                  coords = FALSE)

        init_plot <- ggplot2::ggplot(tmp_data,
                                     aes(x = time,
                                         y = amplitude)) +
          geom_line() +
          facet_grid(electrode ~ epoch,
                     #scales = "free_y",
                     switch = "y") +
          coord_cartesian(ylim = c(-input$uV_scale, input$uV_scale),
                          expand = FALSE) +
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
          geom_vline(xintercept = max(unique(tmp_data$time))) +
          geom_vline(xintercept = 0,
                     linetype = "longdash",
                     alpha = 0.2) +
          labs(x = "Time (s)")
        init_plot
      },
      height = 2500),
      input$dc_offset,
      tmp_data(),
      input$uV_scale
    )

    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })
    session$onSessionEnded(shiny::stopApp)
  }
  shiny::runGadget(ui,
                   server,
                   viewer = shiny::paneViewer(minHeight = 600))
}

