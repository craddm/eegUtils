#' Plot event-related potentials using a scalp based layout
#'
#' Creates a `ggplot2` figure showing the ERP at each electrode, with each
#' electrode placed according to its location on the scalp.
#'
#' @param data An `eeg_epochs` or `eeg_evoked` object
#' @param ... Various arguments passed to specific functions
#'
#' @examples
#' erp_scalp(demo_epochs, montage = "biosemi64alpha")
#'
#' @author Matti Vuorre, \email{mv2521@@columbia.edu}
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @importFrom purrr map
#' @importFrom dplyr group_by summarise mutate select
#' @import ggplot2
#' @import tidyr
#' @family scalp-based maps
#' @seealso [interactive_scalp()] for interactive plots of ERPs in a scalp-based
#'   layout.
#' @return Returns a `ggplot2` plot object.
#' @export

erp_scalp <- function(data,
                      ...) {
  UseMethod("erp_scalp", data)
}

#' @export
erp_scalp.eeg_ICA <- function(data,
                              ...) {
  stop("eeg_ICA objects cannot be used to create scalp plots, since the components do not correspond to single electrodes")
}

#' @param electrode Column name containing electrode names in data. Defaults to
#'   "electrode".
#' @param amplitude Column name containing amplitudes in data. Defaults to
#'   "amplitude".
#' @param time Column name containing time in data. Defaults to "time".
#' @param colour Variable to colour lines by. If no variable is passed, only one
#'   line is drawn for each electrode.
#' @param color Alias for `colour`.
#' @param size Size of the line(s) for the ERPs. Deprecated
#' @param linewidth Size of the line(s) for the ERPs.
#' @param show_guide Should a guide showing the scale of the ERP plots be shown.
#'   Defaults to TRUE.
#' @param baseline Character vector of times to subtract for baseline correct.
#' @param chan_info Pass channel information in the standard `eegUtils` format
#'   directly.
#' @param montage Name of an existing montage set. Defaults to NULL, which will
#'   attempt to us locations from the data. If none are found, it will attempt
#'   to use standard 10-05 locations.
#' @describeIn erp_scalp Plot an ERP scalp map
#' @export
erp_scalp.default <- function(data,
                              electrode = "electrode",
                              amplitude = "amplitude",
                              time = "time",
                              color = NULL,
                              colour = NULL,
                              size = .8,
                              linewidth = .8,
                              baseline = NULL,
                              show_guide = TRUE,
                              chan_info = NULL,
                              montage = NULL,
                              ...) {

  if (is.eeg_epochs(data) && is.null(montage)) {
    chan_info <- channels(data)
    data <- eeg_average(data)
    data <- as.data.frame(data,
                          long = TRUE)
  } else if (is.eeg_epochs(data)) {
    data <- eeg_average(data)
    data <- as.data.frame(data,
                          long = TRUE)
  }

  color <- rlang::enquo(color)
  colour <- rlang::enquo(colour)

  if (rlang::quo_is_null(colour) && !rlang::quo_is_null(color)) {
    colour <- color
  } else if (!rlang::quo_is_null(color) && !rlang::quo_is_null(colour)) {
    message("Please only supply `colour` or `color`, not both. `colour` will be used.")
  }

  if (!rlang::quo_is_symbol(colour) && !rlang::quo_is_null(colour)) {
    colour <- rlang::ensym(colour)
    colour <- rlang::enquo(colour)
  }

  if (rlang::quo_is_null(colour)) {
    data <- dplyr::group_by(data,
                            electrode,
                            time)
    data <- dplyr::summarise(data,
                             amplitude = mean(amplitude))
  } else {
    data <- dplyr::group_by(data,
                            electrode,
                            time,
                            !!colour)
    data <- dplyr::summarise(data,
                             amplitude = mean(amplitude))
  }


  data <- as.data.frame(data)
  # Data maxima for plot limits
  maxAmp <- max(data[, amplitude])
  minAmp <- min(data[, amplitude])
  maxTime <- max(data[, time])
  minTime <- min(data[, time])

  if (!is.null(baseline)) {
    data <- rm_baseline(data, baseline)
  }

  plotfun <- function(x) {
    elec <- unique(x$electrodefacet)
    plot <-
      ggplot2::ggplot(data = x,
                      aes(time,
                          amplitude)) +
      labs(subtitle = elec)

    plot <-
      plot +
      coord_cartesian(ylim = c(minAmp, maxAmp),
                      xlim = c(minTime, maxTime))
    plot <-
      plot +
      geom_vline(xintercept = 0,
                 linewidth = .3)
    plot <-
      plot +
      geom_hline(yintercept = 0,
                 linewidth = .3)
    plot <-
      plot +
      theme_void()
    plot <-
      plot +
      theme(strip.text = element_text(size = 8),
            title = element_text(size = 8))

    if (rlang::quo_is_null(colour)) {
      plot <-
        plot +
        geom_line(linewidth = linewidth)
    } else {
      plot <-
        plot +
        scale_color_brewer(palette = "Set1")
      plot <-
        plot +
        geom_line(aes(colour = {{ colour }}),
                  show.legend = FALSE,
                  linewidth = linewidth)
    }
    plot
  }

  data$electrodefacet <- data[, electrode]
  data <- tidyr::nest(data,
                      data = !electrode)
  data <- dplyr::mutate(data,
                        plot = purrr::map(data,
                                          plotfun))
  data <- dplyr::select(data,
                        -data)

  # Get default electrode locations from pkg internal data
  if (is.null(chan_info)) {
    data <- electrode_locations(data,
                                drop = TRUE,
                                montage = montage)
  } else {
    data <- dplyr::left_join(data,
                             chan_info,
                             by = "electrode")
    data <- dplyr::filter(data, !(is.na(x) | is.na(y)))
  }

  max_x <- max(abs(min(data$x, na.rm = TRUE)),
               abs(max(data$x, na.rm = TRUE)))
  max_y <- max(abs(min(data$y, na.rm = TRUE)),
               abs(max(data$y, na.rm = TRUE)))

  plot_area <- max_x * max_y

  panel_size <- floor(sqrt(plot_area / nrow(data)))

  coords_space <-
    data.frame(x = c(-(max_x + panel_size),
                     max_x + panel_size),
               y = c(-(max_y + panel_size),
                     max_y + panel_size))

  p <-
    ggplot(coords_space,
           aes(x, y)) +
    geom_blank() +
    theme_void() +
    theme(plot.margin = unit(c(10, 10, 10, 10), "pt"))

  guide <-
    ggplot(data,
           aes(x = time,
               y = amplitude)) +
    coord_cartesian(ylim = c(minAmp,
                             maxAmp),
                    xlim = c(minTime,
                             maxTime)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = c(round(minAmp, 1), 0, round(maxAmp, 1))) +
    geom_vline(xintercept = 0,
               linewidth = .4) +
    geom_hline(yintercept = 0,
               linewidth = .4) +
    # Remove the labels, because they make the guide appear TINY
    labs(y = "", x = "") +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = .3),
      axis.title = element_blank(),
      panel.background = element_blank()
    )

  plot_area <- max_x * max_y

  panel_size <- floor(sqrt(plot_area / nrow(data)))

  if (show_guide) {
    p <- p +
      annotation_custom(
        grob = ggplotGrob(guide),
        xmin = min(data$x) - panel_size,
        xmax = min(data$x) + panel_size,
        ymin = min(data$y) - panel_size,
        ymax = min(data$y) + panel_size
      )
  }

  for (i in seq_len(nrow(data))) {
    p <-
      p +
      annotation_custom(
        grob = ggplotGrob(data$plot[[i]]),
        xmin = data$x[i] - panel_size,
        xmax = data$x[i] + panel_size,
        ymin = data$y[i] - panel_size,
        ymax = data$y[i] + panel_size
      )
  }
  p
}

#' Interactive scalp maps
#'
#' Launches a Shiny Gadget for an interactive version of `erp_scalp`, allowing
#' clicking of individual electrodes to plot them in a separate panel. In that
#' panel, they can be averaged over or plotted as individual electrodes.
#'
#' @param data An EEG dataset.
#' @param ... Additional arguments
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @import shiny
#' @import miniUI
#' @family scalp-based maps
#' @seealso [erp_scalp()] for non-interactive plots of ERPs in a scalp-based
#'   layout.
#' @export
interactive_scalp <- function(data,
                              ...) {
  UseMethod("interactive_scalp", data)
}

#' @export
interactive_scalp.default <- function(data,
                                      ...) {
  message(paste("Not implemented for objects of class", class(data)))
}

#' @param colour Variable to colour lines by. If no variable is passed, only one
#'   line is drawn for each electrode.
#' @param baseline Character vector of times to subtract for baseline
#'   correction.
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#' @describeIn interactive_scalp Method for `eeg_epochs` objects
#' @export
interactive_scalp.eeg_epochs <- function(data,
                                         colour = NULL,
                                         baseline = NULL,
                                         montage = NULL,
                                         ...) {

  data <- eeg_average(data, verbose = FALSE)

  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline)
  }

  colour <- rlang::enquo(colour)

  if (!rlang::quo_is_symbol(colour) && !rlang::quo_is_null(colour)) {
    colour <- rlang::ensym(colour)
    colour <- rlang::enquo(colour)
  }

  ui <- miniPage(gadgetTitleBar("Scalp ERPs"),
    miniTabstripPanel(
      miniTabPanel(
        title = "Whole scalp",
        icon = icon("fa-circle", class = "fas"),
        miniContentPanel(
          fillCol(
            flex = c(7, 1),
            plotOutput("Scalp", height = "100%",
                       click = "click_plot"),
            verbatimTextOutput("click_info")
          )
        ),
        miniButtonBlock(actionButton("reset", "Reset selection"))
      ),
      miniTabPanel(
        "Selected electrodes",
        icon = icon("chart-line"),
        miniContentPanel(plotOutput("Selected", height = "100%")),
        miniButtonBlock(
          actionButton("avg", "Mean"),
          actionButton("single", "Plot individual electrodes")
        )
      )
    )
  )

  server <- function(input,
                     output,
                     session) {

    chan_info <- channels(data)

    tmp_data <- as.data.frame(data,
                              long = TRUE)

    button_reacts <- reactiveValues(sel_elecs = list(),
                                    avg = TRUE)

    output$Scalp <- renderPlot({

      if (rlang::quo_is_null(colour)) {
        erp_scalp(data,
                  montage = montage,
                  chan_info = chan_info)
      } else {
        erp_scalp(
          data,
          colour = !!colour,
          montage = montage,
          chan_info = chan_info
        )
      }
    })

    observeEvent(input$click_plot, {
      electrode_choice <- nearPoints(
        tmp_data,
        input$click_plot,
        "x",
        "y",
        threshold = 45,
        maxpoints = 1
      )

      if (nrow(electrode_choice) > 0) {
        if (electrode_choice$electrode %in% button_reacts$sel_elecs) {
          button_reacts$sel_elecs <-
            button_reacts$sel_elecs[-which(button_reacts$sel_elecs == electrode_choice$electrode)]
        } else {
          button_reacts$sel_elecs <- c(button_reacts$sel_elecs,
                                       electrode_choice$electrode)
        }
      }

      output$click_info <- shiny::renderPrint({
        cat("Selected:", unlist(button_reacts$sel_elecs))
      })

      # plot selected electrodes when on appropriate tab
      output$Selected <- shiny::renderPlot({
        if (button_reacts$avg) {
          if (rlang::quo_is_null(colour)) {
            plot_timecourse(data,
                            unlist(button_reacts$sel_elecs))
          } else {
            plot_timecourse(data,
                            unlist(button_reacts$sel_elecs),
                            mapping = aes(colour = !!colour))
          }
        } else {
          if (rlang::quo_is_null(colour)) {
            plot_timecourse(data,
                            electrode = unlist(button_reacts$sel_elecs)) +
              aes(colour = electrode)
          } else {
            plot_timecourse(data,
                            electrode = unlist(button_reacts$sel_elecs),
                            mapping = aes(colour = !!colour)) +
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
