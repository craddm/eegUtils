#' Artefact browser
#'
#' An interactive Shiny app that allows exploration of channel and epoch
#' statistics.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @import shiny
#' @import ggplot2
#' @importFrom plotly plot_ly renderPlotly event_data
#' @param data Object to be explored.
#' @return Nothing.
#' @examples
#'   \dontrun{ view_artefacts(demo_epochs) }
#' @export

view_artefacts <- function(data) {

  chan_dat <- channel_stats(data)
  scale_chans <- dplyr::mutate_if(chan_dat,
                                  is.numeric,
                                  ~(. - mean(.)) / sd(.))
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
    navbarPage("Artefact checks",
               inverse = TRUE,
               tabPanel("Channel stats",
                        sidebarLayout(
                          sidebarPanel(selectInput("chan_meas",
                                                   "Display measures",
                                                   choices = c("means",
                                                               "sds",
                                                               "variance",
                                                               "kurtosis",
                                                               "minmax")
                                                   ),
                                       checkboxInput("std_meas",
                                                     "Standardize?"),
                                       width = 3),
                          mainPanel(
                            plotly::plotlyOutput("chan_plot"),
                            #channelPlotly(chan_dat),
                            fluidRow(
                              column(6,
                                     plotly::plotlyOutput("erpplot")),
                              column(6,
                                     plotly::plotlyOutput("erpimage"))),
                            ),
                          )),
               tabPanel("Epoch stats",
                        epochPlotly(epoch_dat))
    )

  server <- function(input, output) {

    plot_data <- reactive({
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
          labs(y = "amplitude")
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


epochPlotly <- function(id,
                        label = "Plotly channels") {

  ns <- shiny::NS(id)

  tagList(
    fluidRow(
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


