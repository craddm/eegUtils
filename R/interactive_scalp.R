#' Interactive scalp maps
#'
#' Launchs a Shiny Gadget for an interactive version of erp_scalp, allowing clicking of individual electrodes to plot them separately.
#'
#' @param data An EEG dataset.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#'
#' @import shiny
#' @import miniUI
#' @export

interactive_scalp <- function(data) {

  ui <- miniPage(
    gadgetTitleBar("Scalp ERPs"),
    miniTabstripPanel(
      miniTabPanel("Whole scalp", icon = icon("sliders"),
        miniContentPanel(
          fillCol(
            flex = c(4,1),
            plotOutput("Scalp", height = "100%",
                       click = "click_plot"),
            verbatimTextOutput("click_info")
          )
        )
      ),
      miniTabPanel("Selected electrodes", icon = icon("line-chart"),
        miniContentPanel(
          plotOutput("Selected", height = "100%")
          )
        )
      )
    )

  server <- function(input, output, session) {

    data <- electrode_locations(data)
    elec_list <- reactiveValues(sel_elecs = NULL)

    output$Scalp <- renderPlot({
      erp_scalp(data)
      })

    observeEvent(input$click_plot, {
      tmp <- nearPoints(data, input$click_plot, "x", "y", threshold = 45, maxpoints = 1)

      if (tmp$electrode %in% elec_list) {
        elec_list <- elec_list[-which(elec_list == tmp$electrode)]
      } else{
        elec_list <- c(elec_list,tmp$electrode)
      }

      output$click_info <- renderPrint({
        paste("Selected:", elec_list)
        })

      output$Selected <- renderPlot({
        plot_timecourse(data[data$electrode == tmp$electrode, ])
      })
    })

    observeEvent(input$done, {
      stopApp()
    })

  }

  runGadget(ui, server)
}
