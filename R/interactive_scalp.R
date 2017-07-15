#' Interactive scalp maps
#'
#' Launchs a Shiny Gadget for an interactive version of erp_scalp, allowing clicking of individual electrodes to plot them separately.
#'
#' @param data An EEG dataset.
#' @param colour Variable to color lines by. If no variable is passed, only
#' one line is drawn for each electrode.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#'
#' @import shiny
#' @import miniUI
#' @export

interactive_scalp <- function(data, colour = NULL) {

  ui <- miniPage(
    gadgetTitleBar("Scalp ERPs"),
    miniTabstripPanel(
      miniTabPanel("Whole scalp", icon = icon("sliders"),
        miniContentPanel(
          fillCol(
            flex = c(5,1),
            plotOutput("Scalp", height = "100%",
                       click = "click_plot"),
            verbatimTextOutput("click_info")
          )
        ),
        miniButtonBlock(
          actionButton("reset", "Reset selection")
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

    elec_list <- reactiveValues(sel_elecs = list())

    output$Scalp <- renderPlot({
      if (is.null(colour)) {
        erp_scalp(data)
      } else {
        erp_scalp(data, color = as.name(colour))
        }
      })

    observeEvent(input$click_plot, {

      tmp <- nearPoints(data, input$click_plot, "x", "y", threshold = 45, maxpoints = 1)

      if (nrow(tmp) > 0) {
        if (tmp$electrode %in% elec_list$sel_elecs) {
         elec_list$sel_elecs <- elec_list$sel_elecs[-which(elec_list$sel_elecs == tmp$electrode)]
       } else {
         elec_list$sel_elecs <- c(elec_list$sel_elecs, tmp$electrode)
       }
      }

      output$click_info <- renderPrint({
         cat("Selected:", unlist(elec_list$sel_elecs))
        })

      # plot selected electrodes when on appropriate tab
      output$Selected <- renderPlot({
        if (is.null(colour)) {
          plot_timecourse(data[data$electrode %in% elec_list$sel_elecs, ])
        } else{
          plot_timecourse(data[data$electrode %in% elec_list$sel_elecs, ], colour = as.name(colour))
          }
      })
    })

    #Close app gracefully if done button clicked

    observeEvent(input$done, {
      stopApp()
    })

    # Check if reset button clicked on Whole scalp page
    observeEvent(input$reset,{
      elec_list$sel_elecs <- NULL
    })

  }

  runGadget(ui, server)
}
