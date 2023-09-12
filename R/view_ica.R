#' EEG decomposition viewer
#'
#' A Shiny viewer for Independent Component Analysis or Spatio-spectral
#' Decomposition/RESS components that provides an interface for looking at
#' topographies, timecourses, and power spectral densities of all or individual
#' components. Can be used to select and reject artefactual components.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data An `eeg_ICA` object
#' @return A list consisting (optionally) of
#' * A character vector of components marked for rejection
#' * A character vector of components marked to be kept
#' * An `eeg_epochs` object reconstructed from the `eeg_ICA` object, with
#'   components marked for rejection removed.
#' @export

view_ica <- function(data) {

  if (!is.eeg_ICA(data)) {
    stop("This function requires an eeg_ICA object.")
  }

  psd_ica <- compute_psd(data,
                         verbose = FALSE,
                         keep_trials = FALSE)

  psd_ica <- tidyr::pivot_longer(
    psd_ica,
    cols = channel_names(data),
    names_to = "component",
    values_to = "power"
  )

  psd_ica$power <- 10 * log10(psd_ica$power)

  ica_erps <-
    as.data.frame(eeg_average(data,
                              cols = "participant_id"),
                  long = TRUE,
                  coords = FALSE)
  ica_butter <- plot_butterfly(ica_erps)

  ica_topoplots <-
    topoplot(data,
             seq_along(names(data$signals)),
             grid_res = 67,
             chan_marker = "none",
             limits = c(-3, 3),
             verbose = FALSE)

  ui <-
    navbarPage(
      title = "EEG decomposition viewer",
      id = "main_page",
      inverse = TRUE,
      collapsible = TRUE,
      tabPanel(
        "Topographies", plotOutput("comp_topos", dblclick = "topo_click")
      ),
      tabPanel(
        "Timecourses",
        plotOutput("ica_butters",
                   hover = "butter_click",
                   brush = brushOpts(id = "butter_brush",
                                     resetOnNew = TRUE),
                   dblclick = "butter_dbl"),
        tableOutput("info"),
        tags$p(span("Hover over lines to see component details.")),
        tags$p("To zoom, drag-click to highlight where you want to zoom, then double-click to zoom. Double-click again to zoom back out.")
      ),
      tabPanel(
        "PSDs",
        plotOutput(
          "ica_psd",
          hover = "psd_click",
          brush = brushOpts(id = "psd_brush",
                            resetOnNew = TRUE),
          dblclick = "psd_dbl"
        ),
        tableOutput("psd_info"),
        tags$p(span("Hover over lines to see component details.")),
        tags$p("To zoom, drag-click to highlight where you want to zoom, then double-click to zoom. Double-click again to zoom back out.")
      ),
      tabPanel("Individual",
        sidebarLayout(
          sidebarPanel(
                       selectInput(
                         "comp_no", "Component:",
                         channel_names(data)
                       ),
                       shiny::radioButtons("reject_comps",
                                           label = NULL,
                                           choices = c("Keep", "Reject"),
                                           inline = TRUE),
                       width = 3),
          mainPanel(
                    fluidRow(
                      column(plotOutput("indiv_topo"), width = 6),
                      column(plotOutput("indiv_erpim"), width = 6)
                    ),
                    fluidRow(
                      column(plotOutput("indiv_psd"), width = 6),
                      column(plotOutput("indiv_tc"), width = 6)
                    ),
                    width = 9)
        )
      ),
      tabPanel("Output",
               tableOutput("reject_table"),
               checkboxGroupInput("output_choices",
                 label = "Output to return",
                 choices = list(
                                "Components to reject" = "reject",
                                "Components to keep" = "keep",
                                "Reconstructed data" = "data")
               ),
               actionButton("done",
                            "Press to close app and return to console"))

    )

  server <- function(input,
                     output,
                     session) {

    comp_status <- reactiveValues()
    ranges <- reactiveValues(x = NULL,
                             y = NULL)
    b_ranges <- reactiveValues(x = NULL,
                               y = NULL)

    output$comp_topos <-
      shiny::renderPlot(
        ica_topoplots,
        height = function() {
          .75 * session$clientData$output_comp_topos_width
        }
      )

    output$ica_butters <-
      renderPlot(
        ica_butter +
          coord_cartesian(xlim = b_ranges$x,
                          ylim = b_ranges$y)
      )

    output$ica_psd <-
      renderPlot({
        ggplot(psd_ica,
               aes(x = frequency,
                   y = power,
                   colour = component)) +
          geom_line() +
          theme_bw() +
          coord_cartesian(xlim = ranges$x,
                          ylim = ranges$y,
                          expand = FALSE)
      })

    output$info <- renderTable({
      as.data.frame(
        shiny::nearPoints(ica_erps,
                          input$butter_click,
                          threshold = 20,
                          maxpoints = 1,
                          xvar = "time",
                          yvar = "amplitude")
      )
    })

    output$psd_info <- renderTable({
      as.data.frame(
        shiny::nearPoints(psd_ica,
                          input$psd_click,
                          threshold = 20,
                          maxpoints = 1)
      )
    })

    output$indiv_topo <- renderCachedPlot({
      topoplot(data,
               input$comp_no,
               verbose = FALSE)
    }, cacheKeyExpr = {
      input$comp_no
    })

    output$indiv_erpim <- renderCachedPlot({
      erp_image(data,
                input$comp_no)
    }, cacheKeyExpr = {
      input$comp_no
    })

    output$indiv_tc <- renderCachedPlot({
      plot_timecourse(data,
                      input$comp_no)
    },
    cacheKeyExpr = {
      input$comp_no
    })

    output$indiv_psd <- renderCachedPlot({
      tmp_psd <-
        compute_psd(
          select(data,
                 input$comp_no),
          n_fft = data$srate,
          noverlap = 0,
          verbose = FALSE
        )

      tmp_psd <- dplyr::rename(tmp_psd,
                               power = 2)
      tmp_psd <- dplyr::filter(tmp_psd,
                               frequency >= 3,
                               frequency <= 50)
      ggplot(tmp_psd,
             aes(x = frequency,
                 y = 10 * log10((power)))) +
        stat_summary(geom = "ribbon",
                     fun.data = mean_cl_normal,
                     alpha = 0.5) +
        stat_summary(geom = "line",
                     fun = mean) +
        theme_classic() +
        labs(x = "Frequency (Hz)",
             y = "Power (dB)") +
        coord_cartesian(expand = FALSE)
    },
    cacheKeyExpr = {
      input$comp_no
    },
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

    observeEvent(input$topo_click, {
      selected_topo <- as.data.frame(
        shiny::nearPoints(ica_topoplots$data,
                          input$topo_click,
                          threshold = 20,
                          maxpoints = 1)
      )
      updateNavbarPage(inputId = "main_page",
                       selected = "Individual")
      updateSelectInput(inputId = "comp_no",
                        selected = selected_topo$component)
    })


    observeEvent(input$comp_no, {
      updateRadioButtons(inputId = "reject_comps",
                         choices = c("Keep", "Reject"),
                         selected = comp_status[[input$comp_no]],
                         inline = TRUE)
    })

    observeEvent(input$reject_comps, {
      comp_status[[isolate(input$comp_no)]] <- input$reject_comps
      comp_status
    })

    output$reject_table <- renderTable({
      rejects <- reactiveValuesToList(comp_status)
      rejects <-
        names(rejects)[vapply(rejects,
                              function(x) identical(x, "Reject"),
                              logical(1))]
      data.frame("Rejected" = rejects)
    }
    )

    observeEvent(input$done, {
      outputs <- isolate(input$output_choices)

      if (is.null(outputs)) {
        message("No output requested.")
        shiny::stopApp()
      } else {
        returnValue <- vector(
          "list",
          length(outputs)
        )
        names(returnValue) <- outputs
        rejects <- reactiveValuesToList(isolate(comp_status))
        rejects <-
          names(rejects)[vapply(rejects,
                                function(x) identical(x, "Reject"),
                                logical(1))]
        if ("reject" %in% outputs) {
          returnValue$reject <- rejects
        }

        if ("data" %in% outputs) {
          returnValue$data <-
            apply_ica(data,
                      rejects)
        }

        if ("keep" %in% outputs) {
          returnValue$keep <-
            channel_names(data)[!(channel_names(data) %in% rejects)]
        }
      }
      shiny::stopApp(returnValue)

    })
  }
  shiny::runGadget(ui,
                   server)
}


ica_topos <- function(data) {
  ggplot2::ggplot(get_scalpmap(data,
                               grid_res = 50),
                  aes(
                    x = x,
                    y = y,
                    fill = scale(fill),
                    z = scale(fill)
                  )) +
    geom_raster(interpolate = TRUE) +
    geom_head(data = channels(data),
              mapping = aes(fill = NULL,
                            z = NULL)) +
    geom_contour(
      aes(linetype = stat(level) < 0),
      bins = 6,
      colour = "black",
      size = rel(0.8)
    ) +
    facet_wrap(~component) +
    theme_void() +
    scale_fill_distiller(palette = "RdBu",
                         limits = c(-3, 3),
                         oob = scales::squish) +
    theme(legend.position = "none") +
    coord_fixed()
}
