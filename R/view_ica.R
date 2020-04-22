#' EEG component viewer
#'
#' A Shiny viewer for ICA or SSD/RESS components that provides an interface for
#' looking at topographies, timecourses, and power spectral densities of all or
#' individual components.
#'
#' @author Matt craddock \email{matt@@mattcraddock.com}
#' @param data An \code{eeg_ICA} object
#' @export

view_ica <- function(data) {

  if (!is.eeg_ICA(data)) {
    stop("This function requires an eeg_ICA object.")
  }

  ui <- navbarPage("ICA viewer",
                   tabPanel("Topographies",
                            plotOutput("comp_topos",
                                       click = "topo_click"
                                       )),
                   tabPanel("Timecourses",
                            plotOutput("ica_butters",
                                       hover = "butter_click",
                                       brush = brushOpts(
                                         id = "butter_brush",
                                         resetOnNew = TRUE
                                       ),
                                       dblclick = "butter_dbl"),
                            verbatimTextOutput("info")),
                   tabPanel("PSDs",
                            plotOutput("ica_psd",
                                       hover = "psd_click",
                                       brush = brushOpts(
                                         id = "psd_brush",
                                         resetOnNew = TRUE
                                       ),
                                       dblclick = "psd_dbl"),
                            verbatimTextOutput("psd_info")),
                   tabPanel("Individual",
                            fluidRow(
                              column(
                                width = 6,
                                plotOutput("indiv_topo")
                              ),
                              column(
                                width = 6,
                                plotOutput("indiv_erpim"))
                              ),
                            fluidRow(
                              column(
                                width = 6,
                                plotOutput("indiv_tc")
                              ),
                              column(
                                width = 6,
                                selectInput("comp_no", "Comp:",
                                            channel_names(data))
                              )
                            )
                   )
                  )

  server <- function(input,
                     output,
                     session) {


    ranges <- reactiveValues(x = NULL,
                             y = NULL)

    ica_erps <- as.data.frame(eeg_average(data),
                              long = TRUE)
    ica_erps <- dplyr::group_by(ica_erps,
                                electrode,
                                time)

    ica_erps <- dplyr::summarise(ica_erps,
                                 amplitude = mean(amplitude))

    psd_ica <- compute_psd(data,
                           n_fft = data$srate * 2,
                           verbose = FALSE)

    psd_ica <- tidyr::pivot_longer(psd_ica,
                                   cols = channel_names(data),
                                   names_to = "component",
                                   values_to = "power")

    psd_ica <- dplyr::group_by(psd_ica,
                               component,
                               frequency)
    psd_ica <- dplyr::summarise(psd_ica,
                                power = mean(10 * log10(power)))

    output$comp_topos <- renderPlot(
     ica_topos(data),
     height = function() {
       .8 * session$clientData$output_comp_topos_width
      }
     )

    output$ica_butters <- renderPlot(
      plot_butterfly(data)
    )

    output$ica_psd <- renderPlot(
      ggplot(psd_ica,
             aes(x = frequency,
                 y = power,
                 colour = component)) +
        geom_line() +
        theme_bw() +
        coord_cartesian(xlim = ranges$x,
                        ylim = ranges$y,
                        expand = FALSE)
    )

    output$info <- renderPrint({
      nearPoints(ica_erps,
                 input$butter_click,
                 threshold = 10,
                 maxpoints = 1)
    })

    output$psd_info <- renderPrint({
      nearPoints(psd_ica,
                 input$psd_click,
                 threshold = 20,
                 maxpoints = 1)
    })

    output$indiv_topo <- renderPlot({
      topoplot(data, input$comp_no)
    })

    output$indiv_erpim <- renderPlot({
      erp_image(data,
                input$comp_no,
                clim = c(-3, 3))
    })

    output$indiv_tc <- renderPlot({
      plot_timecourse(data,
                      input$comp_no)
    })

    observeEvent(input$psd_dbl, {
      brush <- input$psd_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })

  }
  shinyApp(ui, server)

}


ica_topos <- function(data) {
  ggplot2::ggplot(get_scalpmap(data,
                               grid_res = 67),
                  aes(x = x,
                      y = y,
                      fill = scale(fill),
                      z = scale(fill))) +
    geom_raster() +
    geom_head(data = channels(data),
              mapping = aes(fill = NULL,
                            z = NULL)) +
    geom_contour(aes(linetype = stat(level) < 0),
                 bins = 6,
                 colour = "black",
                 size = rel(0.8)) +
    facet_wrap(~component) +
    theme_void() +
    scale_fill_distiller(palette = "RdBu",
                        limits = c(-2.5, 2.5),
                         oob = scales::squish) +
    theme(legend.position = "none") +
    coord_fixed()
}
