#' Artefact browser
#'
#' An interactive Shiny app that allows exploration of channel statistics.
#'
#' @import shiny
#' @import shinydashboard
#' @importFrom tidyr gather spread
#' @import ggplot2
#' @importFrom plotly plotly renderPlotly event_data
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
                                 tabName = "by_chans"),
        shinydashboard::menuItem("Epoch stats",
                                 tabName = "by_epochs"),
        shinydashboard::menuItem("Channel plotly",
                                 tabName = "plotly_chans"),
        shinydashboard::menuItem("Epoch plotly",
                                 tabName = "plotly_epochs"))),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "by_chans",
                                h2("Channel statistics"),
                                channelPlots(chan_dat)
                                # verbatimTextOutput("info")
                                ),
        shinydashboard::tabItem(tabName = "by_epochs",
                                h2("Epoch statistics"),
                                epochPlots(epoch_dat)
                                ),
        shinydashboard::tabItem(tabName = "plotly_chans",
                                h2("plotly chans"),
                                channelPlotly(chan_dat),
                                plotlyOutput("erpplot")
                                ),
        shinydashboard::tabItem(tabName = "plotly_epochs",
                                h2("plotly epochs"),
                                epochPlotly(epoch_dat))
                                #verbatimTextOutput("selection"))
      )
    )
  )

  server <- function(input, output) {

    output$chan_means <- renderPlot({
      ggplot(chan_dat,
             aes(x = electrode,
                 y = means)) +
        geom_point() +
        theme_bw()
    })

    output$chan_vars <- renderPlot({
      ggplot(chan_dat,
             aes(x = electrode,
                 y = variance)) +
        geom_point() +
        theme_bw()
    })

    output$chan_kurt <- renderPlot({
      ggplot(chan_dat,
             aes(x = electrode,
                 y = kurtosis)) +
        geom_point() +
        theme_bw()
    })

    output$epo_means <- renderPlot({
      ggplot(epoch_dat,
             aes(x = epoch,
                 y = max)) +
        geom_point() +
        theme_bw()
      })

    output$epo_vars <- renderPlot({
      ggplot(epoch_dat,
             aes(x = epoch,
                 y = variance)) +
        geom_point() +
        theme_bw()
      })

    output$epo_kurt <- renderPlot({
      ggplot(epoch_dat,
             aes(x = epoch,
                 y = kurtosis)) +
        geom_point() +
        theme_bw()
      })

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

    # output$scatterplot <- plotly::renderPlotly({
    #   s <- plotly::event_data("plotly_click", source = "chan_plotly")
    #   if (length(s)) {
    #     vars <- c(s[["x"]], s[["y"]])
    #     d <- setNames(mtcars[vars], c("x", "y"))
    #     yhat <- fitted(lm(y ~ x, data = d))
    #     plot_ly(d, x = ~x) %>%
    #       add_markers(y = ~y) %>%
    #       add_lines(y = ~yhat) %>%
    #       layout(xaxis = list(title = s[["x"]]),
    #              yaxis = list(title = s[["y"]]),
    #              showlegend = FALSE)
    #   } else {
    #     plotly_empty()
    #   }
    # })

#     button_reacts <- reactiveValues(sel_elecs = list())
#
#     observeEvent(input$chan_mean_click, {
#
#       tmp <- nearPoints(chan_dat,
#                         input$chan_mean_click,
#                         threshold = 10,
#                         maxpoints = 1)
#
#       print(button_reacts$sel_elecs)
#       print(tmp$electrode)
#       #print(sel_chans)
#       #print(chan_dat[unlist(button_reacts$sel_elecs), "electrode"])
#       if (nrow(tmp) > 0) {
#         #button_reacts$sel_elecs <-
#          # button_reacts$sel_elecs[-which(sel_chans == tmp$electrode)]
#       #   } else {
#          button_reacts$sel_elecs <- c(button_reacts$sel_elecs, tmp$electrode)
#          print(button_reacts$sel_elecs)
#       #   }
#        }
#     })
#
#     output$info <- renderPrint({
#       #chan_dat[unlist(button_reacts$sel_elecs), ]
#       #cat("Selected:", unlist(button_reacts$sel_elecs))
#     })

    output$chan_plotly <- plotly::renderPlotly({
      plotly::plot_ly(chan_dat,
                      x = ~electrode,
                      y = ~means) %>%
        add_markers()
    })

    output$plotly_cvars <- plotly::renderPlotly({
      plotly::plot_ly(chan_dat,
                      x = ~electrode,
                      y = ~variance) %>%
        add_markers()
    })

    output$plotly_emeans <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~electrode,
                      y = ~max) %>%
        add_markers()
    })

    output$plotly_evars <- plotly::renderPlotly({
      plotly::plot_ly(epoch_dat,
                      x = ~epoch,
                      y = ~variance) %>%
        add_markers()
    })

  }

  shinyApp(ui, server)

}

#' UI function for plotting channel statistics.
#' @noRd
channelPlots <- function(id,
                         label = "Channel Stats") {

  ns <- shiny::NS(id)

  shiny::tagList(
    fluidRow(
      shinydashboard::box(plotOutput("chan_means",
                                     height = 250,
                                     click = "chan_mean_click")),
      shinydashboard::box(plotOutput("chan_vars",
                                     height = 250)),
      shinydashboard::box(plotOutput("chan_kurt",
                                     height = 250))
      )
  )
}

#' UI function for plotting epoch statistics.
#' @noRd
epochPlots <- function(id, label = "Epoch Stats") {

  ns <- shiny::NS(id)

  tagList(
    fluidRow(
      shinydashboard::box(plotOutput("epo_means",
                                     height = 250)),
      shinydashboard::box(plotOutput("epo_vars",
                                     height = 250)),
      shinydashboard::box(plotOutput("epo_kurt",
                                     height = 250))
    )
  )
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
                                               height = 250))
      # shinydashboard::box(plotly::plotlyOutput("chan_kurt",
      #                                height = 250))
    )
  )
}
