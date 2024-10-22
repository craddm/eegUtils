#' Artefact browser
#'
#' An interactive Shiny app that allows exploration of channel and epoch
#' statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @import shiny
#' @import ggplot2
#' @param data Object to be explored.
#' @examples
#'   \dontrun{ view_artefacts(demo_epochs) }
#' @export

view_artefacts <- function(data) {

  chan_dat <- channel_stats(data)
  scale_chans <- dplyr::mutate(chan_dat,
                               dplyr::across(dplyr::where(is.numeric),
                                             ~(. - mean(.)) / sd(.)))
  epoch_dat <- epoch_stats(data)
  epoch_dat <- tidyr::pivot_longer(epoch_dat,
                                   cols = channel_names(data),
                                   names_to = "electrode",
                                   values_to = "value")
  epoch_dat <- tidyr::pivot_wider(epoch_dat,
                                  id_cols = c(epoch, electrode),
                                  names_from = measure,
                                  values_from = value)

  ui <-
    bslib::page_navbar(title = "Artefact checks",
      inverse = TRUE,
      bslib::nav_panel(
        title = "Channel stats",
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            shiny::selectInput("chan_meas",
                               "Display measures",
                               choices = c("means",
                                           "sds",
                                           "variance",
                                           "kurtosis",
                                           "minmax")
                               ),
            shiny::checkboxInput("std_meas",
                                 "Standardize?")),
          plotly::plotlyOutput("chan_plot"),
          plotly::plotlyOutput("erpplot"),
          plotly::plotlyOutput("erpimage")
          )
        ),
      bslib::nav_panel("Epoch stats", epoch_plotly(epoch_dat))
    )

  server <- function(input, output) {

    plot_data <- shiny::reactive({
      if (input$std_meas) {
        scale_chans
      } else {
        chan_dat
      }
    })

    output$chan_plot <- plotly::renderPlotly({
      ylab <- switch(
        input$chan_meas,
        means = "Mean",
        variance = "Var",
        sds = "Std. dev.",
        kurtosis = "Kurtosis",
        minmax = "Max - min"
      )
      plotly::plot_ly(plot_data(),
                      x = ~electrode,
                      y = ~get(input$chan_meas)) %>%
        plotly::add_markers() %>%
        plotly::layout(yaxis = list(title = ylab))
    })

    output$erpplot <- plotly::renderPlotly({
      s <- plotly::event_data("plotly_click")

      if (length(s)) {
        plot_timecourse(data,
                        electrode = s[["x"]]) +
          labs(y = "Amplitude")
      }
    })

    output$erpimage <- plotly::renderPlotly({
      s <- plotly::event_data("plotly_click")

      if (length(s)) {
        erp_image(data,
                  electrode = s[["x"]])
      }
    })

    output$plotly_emeans <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      y = ~electrode,
                      x = ~epoch,
                      z = ~max,
                      type = "heatmap") %>%
        plotly::layout(title = "Max")
    })

    output$plotly_evars <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~epoch,
                      y = ~electrode,
                      z = ~variance,
                      type = "heatmap") %>%
        plotly::layout(title = "Variance")
    })

    output$plotly_kurt <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~epoch,
                      y = ~electrode,
                      z = ~kurtosis,
                      type = "heatmap") %>%
        plotly::layout(title = "Kurtosis")
    })

    output$plotly_minmax <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~epoch,
                      y = ~electrode,
                      z = ~minmax,
                      type = "heatmap") %>%
        plotly::layout(title = "Max-min")
    })

  }

  shiny::runGadget(ui,
                   server)
}


epoch_plotly <- function(id,
                         label = "Plotly channels") {

  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      plotly::plotlyOutput("plotly_emeans",
                           height = 250),
      plotly::plotlyOutput("plotly_evars",
                           height = 250),
      plotly::plotlyOutput("plotly_kurt",
                           height = 250),
      plotly::plotlyOutput("plotly_minmax",
                           height = 250)
    )
  )
}
