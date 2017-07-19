#' Interactive scalp maps
#'
#' Launches a Shiny Gadget for an interactive version of erp_scalp, allowing clicking of individual electrodes to plot them in a separate panel. In that panel, they can be averaged over or plotted as individual electrodes.
#'
#' @param data An EEG dataset.
#' @param colour Variable to colour lines by. If no variable is passed, only
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
      miniTabPanel("Whole scalp", icon = icon("circle"),
        miniContentPanel(
          fillCol(
            flex = c(7,1),
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
          ),
          miniButtonBlock(
            actionButton("avg", "Mean"),
            actionButton("single", "Plot individual electrodes")
          )
        )
      )
    )

  server <- function(input, output, session) {

    data <- electrode_locations(data)

    button_reacts <- reactiveValues(sel_elecs = list(), avg = TRUE)

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
        if (tmp$electrode %in% button_reacts$sel_elecs) {
          button_reacts$sel_elecs <- button_reacts$sel_elecs[-which(button_reacts$sel_elecs == tmp$electrode)]
       } else {
         button_reacts$sel_elecs <- c(button_reacts$sel_elecs, tmp$electrode)
       }
      }

      output$click_info <- renderPrint({
         cat("Selected:", unlist(button_reacts$sel_elecs))
        })

      # plot selected electrodes when on appropriate tab
      output$Selected <- renderPlot({
        if (button_reacts$avg){
          if (is.null(colour)) {
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ])
            } else {
              plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ], colour = as.name(colour))
            }
        } else {
          if (is.null(colour)) {
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ], colour = "electrode")
          } else{
            plot_timecourse(data[data$electrode %in% button_reacts$sel_elecs, ], colour = as.name(colour)) +
              facet_wrap(~electrode)
          }
        }
      })
    })

    #Close app gracefully if done button clicked

    observeEvent(input$done, {
      stopApp()
    })

    # Check if reset button clicked on Whole scalp page
    observeEvent(input$reset, {
      button_reacts$sel_elecs <- NULL
    })

    observeEvent(input$avg, {
      button_reacts$avg <- TRUE
    })

    observeEvent(input$single, {
      button_reacts$avg <- FALSE
    })

  }

  runGadget(ui, server)
}
