#' Topographical Plotting Function for EEG
#'
#' Allows topographical plotting of functional data. Output is a ggplot2 object.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data An EEG dataset. If the input is a data.frame, then it must have
#'   columns x, y, and amplitude at present. x and y are (Cartesian) electrode
#'   co-ordinates), amplitude is amplitude.
#' @param ... Various arguments passed to specific functions
#' @examples
#' topoplot(demo_epochs)
#' topoplot(demo_epochs, time_lim = c(.1, .2))
#' @export
#'
#' @section Notes on usage of Generalized Additive Models for interpolation: The
#'   function fits a GAM using the `gam` function from `mgcv`. Specifically, it fits
#'   a spline using the model function gam(z ~ s(x, y, bs = "ts", k = 40). Using
#'   GAMs for smooths is very much experimental. The surface is produced from
#'   the predictions of the GAM model fitted to the supplied data. Values at
#'   each electrode do not necessarily match actual values in the data:
#'   high-frequency variation will tend to be smoothed out. Thus, the method
#'   should be used with caution.

topoplot <- function(data,
                     ...) {
  UseMethod("topoplot", data)
}

#' @describeIn topoplot Default method for data frames.
#' @export

topoplot.default <- function(data,
                             ...) {
  stop("Not implemented for objects of class ", paste(class(data), collapse = "/"))
}

#' @param time_lim Timepoint(s) to plot. Can be one time or a range to average
#'   over. If none is supplied, the function will average across all timepoints
#'   in the supplied data.
#' @param limits Limits of the fill scale - should be given as a character
#'   vector with two values specifying the start and endpoints e.g. limits =
#'   c(-2,-2). Will ignore anything else. Defaults to the range of the data.
#' @param chanLocs Allows passing of channel locations (see
#'   `electrode_locations`)
#' @param method Interpolation method. "Biharmonic" or "gam". "Biharmonic"
#'   implements the same method used in Matlab's EEGLAB. "gam" fits a
#'   Generalized Additive Model with k = 40 knots. Defaults to biharmonic spline
#'   interpolation.
#' @param r Radius of cartoon head_shape; if not given, defaults to 1.1 * the
#'   maximum y electrode location.
#' @param grid_res Resolution of the interpolated grid. Higher = smoother but
#'   slower.
#' @param palette Defaults to RdBu if none supplied. Can be any from
#'   RColorBrewer or viridis. If an unsupported palette is specified, switches
#'   to Greens.
#' @param interp_limit "skirt" or "head". Defaults to "skirt". "skirt"
#'   interpolates just past the farthest electrode and does not respect the
#'   boundary of the head_shape. "head" interpolates up to the radius of the
#'   plotted head.
#' @param contour Plot contour lines on topography (defaults to TRUE)
#' @param chan_marker Set marker for electrode locations. "point" = point,
#'   "name" = electrode name, "none" = no marker. Defaults to "point".
#' @param quantity Allows plotting of arbitrary quantitative column. Defaults to
#'   amplitude. Use quoted column names. E.g. "p.value", "t_statistic".
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#' @param colourmap Deprecated, use palette instead.
#' @param highlights Electrodes to highlight (in white).
#' @param scaling Scaling multiplication factor for labels and any plot lines.
#'   Defaults to 1.
#' @param groups Column name for groups to retain.
#' @param verbose Warning messages when electrodes do not have locations.
#'   Defaults to TRUE.
#' @param scale_fac Scaling factor for masking ring
#' @import ggplot2
#' @import tidyr
#' @importFrom dplyr group_by summarise ungroup
#' @import scales
#' @importFrom mgcv gam
#' @describeIn topoplot Topographical plotting of data.frames and other non
#'   eeg_data objects.
#' @export

topoplot.data.frame <- function(data,
                                time_lim = NULL,
                                limits = NULL,
                                chanLocs = NULL,
                                method = "Biharmonic",
                                r = NULL,
                                grid_res = 200,
                                palette = "RdBu",
                                interp_limit = "skirt",
                                contour = TRUE,
                                chan_marker = "point",
                                quantity = "amplitude",
                                montage = NULL,
                                colourmap,
                                highlights = NULL,
                                scaling = 1,
                                groups = NULL,
                                verbose = TRUE,
                                scale_fac = 1.02,
                                ...) {

  if (!missing(colourmap)) {
    warning("Argument colourmap is deprecated, please use palette instead.",
            call. = FALSE)
    palette <- colourmap
  }

  # Filter out unwanted timepoints and find nearest time values in the data
  # --------------

   if (!is.null(time_lim)) {
    data <- select_times(data,
                         time_lim)
   }

  # Check for x and y co-ordinates, try to add if not found
  # --------------
  if (all(c("x", "y") %in% names(data))) {
    if (verbose) {
      message("Using electrode locations from data.")
    }
  } else if (!is.null(chanLocs)) {

    if (!all(c("x", "y") %in% names(chanLocs))) {
      stop("No channel locations found in chanLocs.")
    }

    data$electrode <- toupper(data$electrode)
    chanLocs$electrode <- toupper(chanLocs$electrode)
    data <- merge(data, chanLocs)

  } else if ("electrode" %in% names(data)) {
    data <- electrode_locations(data,
                                drop = TRUE,
                                montage = montage)
    if (verbose) {
      message("Attempting to add standard electrode locations...")
    }
  } else {
    stop("Neither electrode locations nor labels found.")
  }

  # Remove channels with no location
  if (any(is.na(data$x))) {
    data <- data[!is.na(data$x), ]
    if (verbose) {
      warning("Removing channels with no location.")
    }
  }

  # Average over all timepoints ----------------------------

  x <- NULL
  y <- NULL
  electrode <- NULL
  if (is.character(groups)) {
    groups <- as.name(groups)
  }

  if (is.character(quantity)) {
    quantity <- as.name(quantity)
  }

  if (!is.null(groups)) {

    data <- dplyr::group_by(data,
                            x,
                            y,
                            electrode,
                            {{groups}})
    data <- dplyr::summarise(data,
                             z = mean({{quantity}},
                                      na.rm = TRUE))
    data <- dplyr::ungroup(data)
    data <- tidyr::nest(data,
                        data = -{{groups}})

  } else {
    if (is.character(quantity)) {
      quantity <- as.name(quantity)
    }
    data <-
      dplyr::summarise(dplyr::group_by(data,
                                       x,
                                       y,
                                       electrode),
                       z = mean({{quantity}},
                                na.rm = TRUE))

    # Cut the data frame down to only the necessary columns, and make sure it has
    # the right names
    data <- data.frame(x = data$x,
                       y = data$y,
                       fill = data$z,
                       electrode = data$electrode)

    data <- dplyr::ungroup(data)
    data <- tidyr::nest(tibble::as_tibble(data),
                        data = everything())
  }

  # Do the interpolation! ------------------------
  data <- tidyr::unnest(data,
                        cols = c(data))

  # Add head and mask to topoplot
  if (is.null(r)) {
    abs_x_max <- max(abs(data$x), na.rm = TRUE)
    abs_y_max <- max(abs(data$y), na.rm = TRUE)
    r <- switch(interp_limit,
                "head" = sqrt(abs_x_max^2 + abs_y_max^2),
                "skirt" = 95) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
  }

  # Create the actual plot -------------------------------
  topo <-
    ggplot2::ggplot(get_scalpmap(data,
                                 interp_limit = interp_limit,
                                 method = method,
                                 grid_res = grid_res,
                                 r = r),
                    aes(x = x,
                        y = y,
                        fill = fill)) +
     geom_raster(interpolate = TRUE,
                 na.rm = TRUE)

  if (contour) {
    topo <-
      topo +
      stat_contour(
        aes(z = fill,
            linetype = after_stat(level) < 0),
        bins = 6,
        colour = "black",
        size = rel(1.1 * scaling),
        show.legend = FALSE
      )
    }

  topo <-
    topo +
    geom_mask(scale_fac = scale_fac,
              size = 5 * scaling) +
    geom_head(r = r,
              size = rel(1.5) * scaling) +
    coord_equal() +
    theme_bw() +
    theme(rect = element_blank(),
          line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    guides(fill = guide_colorbar(title = expression(paste("Amplitude (",
                                                          mu, "V)")),
                                 title.position = "right",
                                 barwidth = rel(1) * scaling,
                                 barheight = rel(6) * scaling,
                                 title.theme = element_text(angle = 270)))

  # Add electrode points or names -------------------
  if (chan_marker == "point") {
    topo <-
      topo +
      ggplot2::annotate("point",
                       x = data$x,
                       y = data$y,
                       colour = "black",
                       size = rel(2 * scaling))
    }  else if (chan_marker == "name") {
      topo <-
        topo +
        ggplot2::annotate("text",
                          x = data$x,
                          y = data$y,
                          label = data$electrode,
                          colour = "black",
                          size = rel(4 * scaling))
    }

  # Highlight specified electrodes
  if (!is.null(highlights)) {
    high_x <- data$x[data$electrode %in% highlights]
    high_y <- data$y[data$electrode %in% highlights]
    topo <- topo +
      annotate("point",
               x = high_x,
               y = high_y,
               colour = "white",
               size = rel(2 * scaling))
  }

  # Set the palette and scale limits ------------------------
  topo <- set_palette(topo,
                      palette,
                      limits)
  topo
}

#' @describeIn topoplot Topographical plotting of `eeg_data` objects.
#' @export

topoplot.eeg_data <- function(data, time_lim = NULL,
                              limits = NULL,
                              chanLocs = NULL,
                              method = "Biharmonic",
                              r = NULL,
                              grid_res = 200,
                              palette = "RdBu",
                              interp_limit = "skirt",
                              contour = TRUE,
                              chan_marker = "point",
                              quantity = "amplitude",
                              montage = NULL,
                              highlights = NULL,
                              scaling = 1,
                              verbose = TRUE,
                              scale_fac = 1.02,
                              ...) {

  if (!is.null(data$chan_info)) {
    chanLocs <- channels(data)
  }

  if (is.null(time_lim)) {
    data <- as.data.frame(data)
    data <- as.data.frame(t(colMeans(data)))
    data <- tidyr::gather(data,
                          electrode,
                          amplitude,
                          -sample,
                          -time)
  } else {
    data <- select_times(data,
                         time_lim)
    data <- as.data.frame(data,
                          long = TRUE)
  }
  topoplot(data,
           time_lim = time_lim,
           limits = limits,
           chanLocs = chanLocs,
           method = method,
           r = r,
           grid_res = grid_res,
           palette = palette,
           interp_limit = interp_limit,
           contour = contour,
           chan_marker = chan_marker,
           quantity = quantity,
           montage = montage,
           highlights = highlights,
           passed = TRUE,
           scaling = scaling,
           verbose = verbose,
           scale_fac = scale_fac)
}


#' @describeIn topoplot Topographical plotting of `eeg_epochs` objects.
#' @export

topoplot.eeg_epochs <- function(data,
                                time_lim = NULL,
                                limits = NULL,
                                chanLocs = NULL,
                                method = "Biharmonic",
                                r = NULL,
                                grid_res = 200,
                                palette = "RdBu",
                                interp_limit = "skirt",
                                contour = TRUE,
                                chan_marker = "point",
                                quantity = "amplitude",
                                montage = NULL,
                                highlights = NULL,
                                scaling = 1,
                                groups = NULL,
                                verbose = TRUE,
                                scale_fac = 1.02,
                                ...) {

  if (!is.null(data$chan_info)) {
    chanLocs <- channels(data)
  }

  # average over epochs first, but preserve conditions
  if ("event_label" %in% names(data$events)) {
    data <- eeg_average(data,
                        cond_label = list_events(data)$event_label)
  } else {
    data <- eeg_average(data)
  }

  data <- as.data.frame(data,
                        long = TRUE)

  topoplot(data,
           time_lim = time_lim,
           limits = limits,
           chanLocs = chanLocs,
           method = method,
           r = r,
           grid_res = grid_res,
           palette = palette,
           interp_limit = interp_limit,
           contour = contour,
           chan_marker = chan_marker,
           quantity = quantity,
           montage = montage,
           highlights = highlights,
           scaling = scaling,
           groups = groups,
           verbose = verbose,
           scale_fac = scale_fac
           )
}


#' @param component Component to plot (numeric)
#' @describeIn topoplot Topographical plot for `eeg_ICA` objects
#' @export
topoplot.eeg_ICA <- function(data,
                             component,
                             time_lim = NULL,
                             limits = NULL,
                             chanLocs = NULL,
                             method = "Biharmonic",
                             r = NULL,
                             grid_res = 200,
                             palette = "RdBu",
                             interp_limit = "skirt",
                             contour = TRUE,
                             chan_marker = "point",
                             quantity = "amplitude",
                             montage = NULL,
                             colourmap,
                             highlights = NULL,
                             scaling = 1,
                             verbose = TRUE,
                             scale_fac = 1.02,
                             ...) {
  if (missing(component)) {
    stop("Component number must be specified for eeg_ICA objects.")
  }

  if (!is.null(time_lim) && verbose) {
    warning("time_lim is ignored for ICA components.")
  }

  chan_info <- data$chan_info
  data <- data.frame(amplitude = data$mixing_matrix[, component],
                      electrode = data$mixing_matrix$electrode)
  topoplot(data,
           chanLocs = chan_info,
           limits = limits,
           interp_limit = interp_limit,
           r = r,
           grid_res = grid_res,
           palette = palette,
           scaling = scaling,
           method = method,
           quantity = quantity,
           montage = montage,
           contour = contour,
           highlights = NULL,
           chan_marker = chan_marker,
           time_lim = NULL,
           verbose = verbose,
           scale_fac = scale_fac)

}

#' @param freq_range Range of frequencies to average over.
#' @describeIn topoplot Topographical plotting of `eeg_tfr` objects.
#' @export

topoplot.eeg_tfr <- function(data,
                             time_lim = NULL,
                             limits = NULL,
                             chanLocs = NULL,
                             method = "Biharmonic",
                             r = NULL,
                             grid_res = 200,
                             palette = "RdBu",
                             interp_limit = "skirt",
                             contour = TRUE,
                             chan_marker = "point",
                             quantity = "power",
                             montage = NULL,
                             highlights = NULL,
                             scaling = 1,
                             freq_range = NULL,
                             verbose = TRUE,
                             ...) {

  if (!is.null(data$chan_info)) {
    chanLocs <- data$chan_info
  }

  if (!is.null(freq_range)) {
    data <- select_freqs(data, freq_range)
  }

  if (data$freq_info$baseline == "none") {
    palette <- "viridis"
  }

  # average over epochs first
  data <- eeg_average(data)

  data <- as.data.frame(data,
                        long = TRUE)

  topoplot(data,
           time_lim = time_lim,
           limits = limits,
           chanLocs = chanLocs,
           method = method,
           r = r,
           grid_res = grid_res,
           palette = palette,
           interp_limit = interp_limit,
           contour = contour,
           chan_marker = chan_marker,
           quantity = quantity,
           montage = montage,
           highlights = highlights,
           scaling = scaling,
           passed = TRUE,
           verbose = verbose)
}

#' Set palette and limits for topoplot
#'
#' @param topo ggplot2 object produced by topoplot command
#' @param palette Requested palette
#' @param limits Limits of colour scale
#' @import ggplot2
#' @keywords internal

set_palette <- function(topo, palette, limits = NULL) {

  if (palette %in% c("magma", "inferno", "plasma",
                  "viridis", "A", "B", "C", "D")) {

    topo <- topo +
      ggplot2::scale_fill_viridis_c(option = palette,
                                  limits = limits,
                                  guide = "colourbar",
                                  oob = scales::squish)
  } else {
    topo <- topo +
      scale_fill_distiller(palette = palette,
                           limits = limits,
                           guide = "colourbar",
                           oob = scales::squish)
  }
  topo
}
