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

  ui <-
    navbarPage("ICA viewer",
               tabPanel("Topographies",
                        plotOutput("comp_topos",
                                   dblclick = "topo_click"
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
                        sidebarLayout(
                          sidebarPanel(selectInput("comp_no",
                                                   "Component:",
                                                   channel_names(data)),
                                       width = 3),
                          mainPanel(
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
                                plotOutput("indiv_psd")
                                ),
                              column(
                                width = 6,
                                plotOutput("indiv_tc")
                                )
                              ),
                            width = 9
                            )
                          )
                        )
               )

  server <- function(input,
                     output,
                     session) {

    ranges <- reactiveValues(x = NULL,
                             y = NULL)
    b_ranges <- reactiveValues(x = NULL,
                               y = NULL)

    ica_erps <- as.data.frame(eeg_average(data),
                              long = TRUE)

    ica_erps <- dplyr::group_by(ica_erps,
                                electrode,
                                time)

    ica_erps <- dplyr::summarise(ica_erps,
                                 amplitude = mean(amplitude))

    psd_ica <- compute_psd(data,
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

    output$comp_topos <- shiny::renderPlot(
     ica_topos(data),
     height = function() {
       .8 * session$clientData$output_comp_topos_width
      }
     )

    output$ica_butters <- renderPlot(
      plot_butterfly(data) +
        coord_cartesian(xlim = b_ranges$x,
                        ylim = b_ranges$y)
    )

    output$ica_psd <- renderPlot({
      ggplot(psd_ica,
             aes(x = frequency,
                 y = power,
                 colour = component)) +
        geom_line() +
        theme_bw() +
        coord_cartesian(xlim = ranges$x,
                        ylim = ranges$y,
                        expand = FALSE)}
    )

    output$info <- renderPrint({
      shiny::nearPoints(ica_erps,
                        input$butter_click,
                        threshold = 10,
                        maxpoints = 1)
    })

    output$psd_info <- renderPrint({
      shiny::nearPoints(psd_ica,
                        input$psd_click,
                        threshold = 20,
                        maxpoints = 1)
    })

    output$indiv_topo <- renderCachedPlot({
      topoplot(data, input$comp_no)
    }, cacheKeyExpr = {input$comp_no})

    output$indiv_erpim <- renderCachedPlot({
      erp_image(data,
                input$comp_no)
    }, cacheKeyExpr = {input$comp_no})

    output$indiv_tc <- renderCachedPlot({
      plot_timecourse(data,
                      input$comp_no)
    },
    cacheKeyExpr = {input$comp_no})

    output$indiv_psd <- renderCachedPlot({
      tmp_psd <-
        compute_psd(select(data,
                           input$comp_no),
                    n_fft = data$srate,
                    noverlap = 0,
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
        labs(x = "Frequency (Hz)",
             y = "Power (dB)") +
        coord_cartesian(expand = FALSE)
    },
    cacheKeyExpr = {input$comp_no},
    res = 96)

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

    observeEvent(input$butter_dbl, {
      brush <- input$butter_brush
      if (!is.null(brush)) {
        b_ranges$x <- c(brush$xmin,
                      brush$xmax)
        b_ranges$y <- c(brush$ymin,
                      brush$ymax)

      } else {
        b_ranges$x <- NULL
        b_ranges$y <- NULL
      }
    })

  }
  shiny::shinyApp(ui,
                  server)

}


ica_topos <- function(data) {
  ggplot2::ggplot(get_scalpmap(data,
                               grid_res = 80),
                  aes(x = x,
                      y = y,
                      fill = scale(fill),
                      z = scale(fill))) +
    geom_raster(interpolate = TRUE) +
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
                        limits = c(-3, 3),
                        oob = scales::squish) +
    theme(legend.position = "none") +
    coord_fixed()
}
