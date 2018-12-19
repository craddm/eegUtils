#' Plot ERPs on the scalp
#'
#' Creates an ERP figure for each electrode and layouts them on the scalp.
#'
#' @param data An EEG dataset.
#' @param electrode Column name containing electrode names in data.
#' Defaults to "electrode".
#' @param amplitude Column name containing amplitudes in data.
#' Defaults to "amplitude".
#' @param time Column name containing time in data. Defaults to "time".
#' @param color Variable to color lines by. If no variable is passed, only
#' one line is drawn for each electrode.
#' @param size Size of line(s).
#' @param show_guide Should a guide be shown.
#' @param baseline Character vector of times to subtract for baseline correct.
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#'
#' @details The function uses default electrode names and locations contained
#' in the package.
#'
#' @examples
#' erp_scalp(demo_epochs, montage = "biosemi64alpha")
#'
#' @author Matti Vuorre, \email{mv2521@@columbia.edu}
#' @importFrom purrr map
#' @importFrom dplyr group_by summarise mutate select
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @return Returns a ggplot2 plot object.
#' @export

erp_scalp <- function(data,
                      electrode = "electrode",
                      amplitude = "amplitude",
                      time = "time",
                      color = NULL,
                      size = .8,
                      baseline = NULL,
                      show_guide = TRUE,
                      montage = NULL) {

  chan_info <- NULL
  if (is.eeg_epochs(data) & is.null(montage)) {
    chan_info <- channels(data)
    data <- as.data.frame(data,
                          long = TRUE)
  } else if (is.eeg_epochs(data)) {
    data <- as.data.frame(data,
                          long = TRUE)
  }

  if (is.null(color)) {
    data <- dplyr::group_by(data,
                            electrode,
                            time)
    data <- dplyr::summarise(data,
                             amplitude = mean(amplitude))
  } else {
    data <- dplyr::group_by_(data,
                             electrode,
                             time,
                             as.name(color))
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
    plot <- ggplot(data = x,
                   aes_(as.name(time),
                        as.name(amplitude)))
    plot <- plot +
      facet_wrap("electrodefacet")
    plot <- plot +
      coord_cartesian(ylim = c(minAmp, maxAmp),
                      xlim = c(minTime, maxTime))
    plot <- plot +
      geom_vline(xintercept = 0,
                 size = .3)
    plot <- plot +
      geom_hline(yintercept = 0,
                 size = .3)
    plot <- plot +
      theme_void()
    plot <- plot +
      theme(strip.text = element_text(size = 8))
    if (!is.null(color)) {
      plot <- plot + scale_color_brewer(palette = "Set1")
      plot <- plot + geom_line(aes_(color = as.name(color)),
                               show.legend = FALSE, size = size)
    } else {
      plot <- plot + geom_line(size = size)
    }
    plot
  }

  data$electrodefacet <- data[, electrode]
  data <- tidyr::nest(tibble::as.tibble(data), -electrode)
  data <- dplyr::mutate(data, plot = map(data, plotfun))
  data <- dplyr::select(data, -data)

  # Get default electrode locations from pkg internal data
  if (is.null(chan_info)){
    data <- electrode_locations(data,
                                drop = TRUE,
                                montage = montage)
  } else {
    data <- dplyr::left_join(data, chan_info)
    data <- filter(data, !(is.na(x) | is.na(y)))
  }

  max_x <- max(abs(min(data$x, na.rm = TRUE)),
               abs(max(data$x, na.rm = TRUE)))
  max_y <- max(abs(min(data$y, na.rm = TRUE)),
               abs(max(data$y, na.rm = TRUE)))

  plot_area <- max_x * max_y

  panel_size <- floor(sqrt(plot_area / nrow(data)))

  coords_space <- data.frame(x = c(-(max_x + panel_size), max_x + panel_size),
                             y = c(-(max_y + panel_size), max_y + panel_size))

  p <- ggplot(coords_space,
              aes(x, y)) +
    geom_blank() +
    theme_void() +
    theme(plot.margin = unit(c(10, 10, 10, 10), "pt"))

  guide <- ggplot(data,
                  aes(x = time,
                      y = amplitude)) +
    coord_cartesian(ylim = c(minAmp,
                             maxAmp),
                    xlim = c(minTime,
                             maxTime)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    geom_vline(xintercept = 0,
               size = .4) +
    geom_hline(yintercept = 0,
               size = .4) +
    labs(y = expression(paste("Amplitude (",
                              mu,
                              "V)")),
         x = "Time (s)") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(size = .3),
          plot.margin = unit(c(8, 8, 8, 8), "pt"))

  plot_area <- max_x * max_y

  panel_size <- floor(sqrt(plot_area / nrow(data)))

  if (show_guide) {
    p <- p +
      annotation_custom(grob = ggplotGrob(guide),
                        xmin = min(data$x) - panel_size,
                        xmax = min(data$x) + panel_size,
                        ymin = min(data$y) - panel_size,
                        ymax = min(data$y) + panel_size)
  }

  for (i in 1:nrow(data)) {
    p <- p +
      annotation_custom(grob = ggplotGrob(data$plot[[i]]),
                        xmin = data$x[i] - panel_size,
                        xmax = data$x[i] + panel_size,
                        ymin = data$y[i] - panel_size,
                        ymax = data$y[i] + panel_size)
  }

  p
}

#' Interactive scalp maps
#'
#' Launches a Shiny Gadget for an interactive version of erp_scalp, allowing
#' clicking of individual electrodes to plot them in a separate panel. In that
#' panel, they can be averaged over or plotted as individual electrodes.
#' Currently buggy with datasets with existing electrode co-ordinates.
#'
#' @param data An EEG dataset.
#' @param colour Variable to colour lines by. If no variable is passed, only one
#'   line is drawn for each electrode.
#' @param baseline Character vector of times to subtract for baseline correction.
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#'
#' @import shiny
#' @import miniUI
#' @export

interactive_scalp <- function(data,
                              colour = NULL,
                              baseline = NULL,
                              montage = NULL) {

  if (is.eeg_data(data)) {
    data <- eeg_average(data)
  }

  if (!is.null(baseline)) {
    data <- rm_baseline(data,
                        time_lim = baseline)
  }

  ui <- miniPage(
    gadgetTitleBar("Scalp ERPs"),
    miniTabstripPanel(
      miniTabPanel("Whole scalp", icon = icon("circle"),
                   miniContentPanel(
                     fillCol(
                       flex = c(7, 1),
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

  server <- function(input,
                     output,
                     session) {

    chan_info <- NULL
    if (is.eeg_evoked(data)) {
      chan_info <- channels(data)
    }
    tmp_data <- as.data.frame(data,
                              long = TRUE)

    if (is.null(chan_info)){
      tmp_data <- electrode_locations(tmp_data,
                                      drop = TRUE,
                                      montage = montage)
    } else {
      tmp_data <- dplyr::left_join(tmp_data, chan_info)
      tmp_data <- filter(tmp_data, !(is.na(x) | is.na(y)))
    }
    button_reacts <- reactiveValues(sel_elecs = list(),
                                    avg = TRUE)

    output$Scalp <- renderPlot({

      if (is.null(colour)) {
        erp_scalp(tmp_data,
                  montage = montage)
      } else {
        erp_scalp(tmp_data,
                  color = as.name(colour),
                  montage = montage)
      }

    })

    observeEvent(input$click_plot, {

      tmp <- nearPoints(tmp_data,
                        input$click_plot,
                        "x", "y",
                        threshold = 45,
                        maxpoints = 1)

      if (nrow(tmp) > 0) {
        if (tmp$electrode %in% button_reacts$sel_elecs) {
          button_reacts$sel_elecs <-
            button_reacts$sel_elecs[-which(button_reacts$sel_elecs == tmp$electrode)]
        } else {
          button_reacts$sel_elecs <- c(button_reacts$sel_elecs, tmp$electrode)
        }
      }

      output$click_info <- renderPrint({
        cat("Selected:", unlist(button_reacts$sel_elecs))
      })

      # plot selected electrodes when on appropriate tab
      output$Selected <- renderPlot({
        if (button_reacts$avg) {

          if (is.null(colour)) {
            plot_timecourse(select_elecs(data,
                                         unlist(button_reacts$sel_elecs)))
          } else {
            plot_timecourse(select_elecs(data,
                                 unlist(button_reacts$sel_elecs)),
                    colour = as.name(colour))
          }

        } else {
          if (is.null(colour)) {
            plot_timecourse(select_elecs(data,
                                 unlist(button_reacts$sel_elecs)),
                    colour = "electrode")
          } else{
            plot_timecourse(select_elecs(data,
                                         unlist(button_reacts$sel_elecs)),
                             colour = as.name(colour)) +
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


