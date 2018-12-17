#' Artefact browser
#'
#' An interactive Shiny app that allows exploration of channel and epoch statistics.
#'
#' @import shiny
#' @import shinydashboard
#' @importFrom tidyr gather spread
#' @import ggplot2
#' @importFrom plotly plot_ly renderPlotly event_data
#' @keywords internal

view_artefacts <- function(data) {

  chan_dat <- channel_stats(data)
  epoch_dat <- epoch_stats(data)
  epoch_dat <- tidyr::gather(epoch_dat,
                             electrode,
                             value,
                             -measure,
                             -epoch)
  epoch_dat <- tidyr::spread(epoch_dat,
                             measure,
                             value)

  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "Artefact viewer"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Channel stats",
                                 tabName = "plotly_chans"),
        shinydashboard::menuItem("Epoch stats",
                                 tabName = "plotly_epochs"))),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "plotly_chans",
                                h2("plotly chans"),
                                channelPlotly(chan_dat),
                                plotly::plotlyOutput("erpplot")
                                ),
        shinydashboard::tabItem(tabName = "plotly_epochs",
                                h2("plotly epochs"),
                                epochPlotly(epoch_dat))
                                #verbatimTextOutput("selection"))
      )
    )
  )

  server <- function(input, output) {


    output$selection <- renderPrint({
      s <- plotly::event_data("plotly_click")
      if (length(s) == 0) {
        "Click on a cell in the heatmap to display a scatterplot"
      } else {
        cat("You selected: \n\n")
        as.list(s)
      }
    })

    output$erpplot <- plotly::renderPlotly({
      s <- plotly::event_data("plotly_click")

      print(s)
      if (length(s)) {
        plot_timecourse(data,
                        electrode = s[["x"]]) +
          labs(y = "amplitude")
      }
    })

    output$chan_plotly <- plotly::renderPlotly({
      plotly::plot_ly(chan_dat,
                      x = ~electrode,
                      y = ~means) %>%
        plotly::add_markers()
    })

    output$plotly_cvars <- plotly::renderPlotly({
      plotly::plot_ly(chan_dat,
                      x = ~electrode,
                      y = ~variance) %>%
        plotly::add_markers()
    })

    output$plotly_emeans <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      y = ~electrode,
                      x = ~epoch,
                      z = ~max,
                      type = "heatmap")
    })

    output$plotly_evars <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~epoch,
                      y = ~electrode,
                      z = ~variance,
                      type = "heatmap")
    })

    output$plotly_kurt <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~epoch,
                      y = ~electrode,
                      z = ~kurtosis,
                      type = "heatmap")
    })

  }

  shinyApp(ui, server)

}


channelPlotly <- function(id,
                          label = "Plotly channels") {

  ns <- shiny::NS(id)

  tagList(
    fluidRow(
      shinydashboard::box(plotly::plotlyOutput("chan_plotly",
                                     height = 250)),
      shinydashboard::box(plotly::plotlyOutput("plotly_cvars",
                                               height = 250))
      # shinydashboard::box(plotly::plotlyOutput("chan_kurt",
      #                                height = 250))
    )
  )
}


epochPlotly <- function(id,
                        label = "Plotly channels") {

  ns <- shiny::NS(id)

  tagList(
    fluidRow(
      shinydashboard::box(plotly::plotlyOutput("plotly_emeans",
                                               height = 250)),
      shinydashboard::box(plotly::plotlyOutput("plotly_evars",
                                               height = 250)),
      shinydashboard::box(plotly::plotlyOutput("plotly_kurt",
                                    height = 250))
    )
  )
}
