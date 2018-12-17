#' Topographical Plotting Function for EEG
#'
#' Allows topographical plotting of functional data. Output is a ggplot2 object.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data An EEG dataset. If the input is a data.frame, then it must have
#'   columns x, y, and amplitude at present. x and y are (Cartesian) electrode
#'   co-ordinates), amplitude is amplitude.
#' @param ... Various arguments passed to specific functions
#' @export
#'
#' @section Notes on usage of Generalized Additive Models for interpolation: The
#'   function fits a GAM using the gam function from mgcv. Specifically, it fits
#'   a spline using the model function gam(z ~ s(x, y, bs = "ts", k = 40). Using
#'   GAMs for smooths is very much experimental. The surface is produced from
#'   the predictions of the GAM model fitted to the supplied data. Values at
#'   each electrode do not necessarily match actual values in the data:
#'   high-frequency variation will tend to be smoothed out. Thus, the method
#'   should be used with caution.

topoplot <- function(data, ...) {
  UseMethod("topoplot", data)
}

#' @describeIn topoplot Default method for data frames.
#' @export

topoplot.default <- function(data, ...) {
  stop("This function requires a data frame or an eeg_data/eeg_epochs object")
}

#' @param time_lim Timepoint(s) to plot. Can be one time or a range to average
#'   over. If none is supplied, the function will average across all timepoints
#'   in the supplied data.
#' @param limits Limits of the fill scale - should be given as a character vector
#'   with two values specifying the start and endpoints e.g. limits = c(-2,-2).
#'   Will ignore anything else. Defaults to the range of the data.
#' @param chanLocs Allows passing of channel locations (see \code{electrode_locations})
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
#' @param highlights Electrodes to highlight (in white)
#' @param scaling Scaling multiplication factor for labels and any plot lines. Defaults to 1.
#' @param groups Column name for groups to retain.
#' @import ggplot2
#' @import tidyr
#' @importFrom dplyr group_by mutate summarise ungroup
#' @importFrom rlang parse_quo current_env
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
                                grid_res = 67,
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
                                ...) {

  if (!missing(colourmap)) {
    warning("Argument colourmap is deprecated, please use palette instead.",
            call. = FALSE)
    palette <- colourmap
  }

  # Filter out unwanted timepoints, and find nearest time values in the data
  # --------------

   if (!is.null(time_lim)) {
  #   if (length(time_lim) == 1) {
  #     time_lim <- c(time_lim, time_lim)
  #   }
    data <- select_times(data, time_lim)
   }

  # Check for x and y co-ordinates, try to add if not found --------------
  if (!is.null(chanLocs)) {

    if (!all(c("x", "y") %in% names(chanLocs))) {
      stop("No channel locations found in chanLocs.")
    }

    data$electrode <- toupper(data$electrode)
    chanLocs$electrode <- toupper(chanLocs$electrode)
    data <- merge(data, chanLocs)

    # Remove channels with no location
    if (any(is.na(data$x))) {
      data <- data[!is.na(data$x), ]
      warning("Removing channels with no location.")
    }

  } else if (all(c("x", "y") %in% names(data))) {
    message("Using electrode locations from data.")
  } else if ("electrode" %in% names(data)) {
    data <- electrode_locations(data,
                                drop = TRUE,
                                montage = montage)
    message("Attempting to add standard electrode locations...")
  } else {
    stop("Neither electrode locations nor labels found.")
  }


  # Average over all timepoints ----------------------------

  x <- NULL
  y <- NULL
  electrode <- NULL
  if (!is.null(groups)) {
    groups <- rlang::parse_quo(groups,
                               rlang::current_env())
    data <- dplyr::group_by(data,
                            x,
                            y,
                            electrode,
                            !!groups)
    data <- dplyr::summarise(data,
                             z = mean(!!rlang::parse_quo(quantity,
                                                         rlang::current_env()),
                                      na.rm = TRUE))
    data <- dplyr::ungroup(data)
    data <- tidyr::nest(data, -!!groups)
    # Rescale electrode co-ordinates to be from -1 to 1 for plotting
    # Selects largest absolute value from x or y
    max_dim <- max(abs(data$data[[1]]$x),
                   abs(data$data[[1]]$y))
    scaled_x <- data$data[[1]]$x / max_dim
    scaled_y <- data$data[[1]]$y / max_dim
  } else {
    data <- dplyr::summarise(dplyr::group_by(data,
                                             x,
                                             y,
                                             electrode),
                             z = mean(!!rlang::parse_quo(quantity,
                                                         rlang::current_env()),
                                      na.rm = TRUE))
    # Cut the data frame down to only the necessary columns, and make sure it has
    # the right names
    data <- data.frame(x = data$x,
                       y = data$y,
                       z = data$z,
                       electrode = data$electrode)
    # Rescale electrode co-ordinates to be from -1 to 1 for plotting
    # Selects largest absolute value from x or y
    max_dim <- max(abs(data$x),
                   abs(data$y))
    scaled_x <- data$x / max_dim
    scaled_y <- data$y / max_dim
    data <- dplyr::ungroup(data)
    data <- tidyr::nest(tibble::as_tibble(data))
  }

  # Do the interpolation! ------------------------

  switch(method,
         Biharmonic = {
           out_df <-
             dplyr::mutate(data,
                           topos = map(data,
                                       ~biharm_topo(.x,
                                                    grid_res = grid_res,
                                                    scaled_x = scaled_x,
                                                    scaled_y = scaled_y)))
         },
         gam = {
           out_df <- dplyr::mutate(data,
                                 topos = map(data,
                                             ~gam_topo(.x,
                                                       grid_res = grid_res,
                                                       scaled_x = scaled_x,
                                                       scaled_y = scaled_y)))
         })

  out_df <- tidyr::unnest(out_df,
                          topos)

  data <- tidyr::unnest(data)

  # Create the head_shape -----------------

  #set radius as max of y (i.e. furthest forward electrode's y position). Add a
  #little to push the circle out a bit more. Consider making max of abs(scaled_y)

  if (is.null(r)) {
    r <- max(scaled_y) * 1.1
  }

  circ_rads <- seq(0, 2 * pi, length.out = 101)

  head_df <- create_head(r, circ_rads)
  head_shape <- head_df$head_shape
  nose <- head_df$nose
  ears <- head_df$ears

  # Check if should interp/extrap beyond head_shape, and set up ring to mask
  # edges for smoothness

  out_df <- round_topo(out_df,
                       r = r,
                       interp_limit = interp_limit,
                       circ_rads = circ_rads)
  mask_ring <- out_df$mask_ring
  out_df <- out_df$out_df

  # Create the actual plot -------------------------------

  topo <- ggplot2::ggplot(out_df,
                          aes(x,
                              y,
                              fill = amplitude)) +
    geom_raster(interpolate = TRUE, na.rm = TRUE)

  if (contour) {
    topo <- topo +
      stat_contour(
        aes(z = amplitude,
            linetype = ..level.. < 0),
        bins = 6,
        colour = "black",
        size = rel(1.1 * scaling),
        show.legend = FALSE
      )
    }

  # Add head and mask to topoplot
  topo <- topo +
     annotate("path",
              x = mask_ring$x,
              y = mask_ring$y,
              colour = "white",
              size = rel(6.5)) +
    annotate("path",
             x = head_shape$x,
             y = head_shape$y,
             size = rel(1.5 * scaling)) +
    annotate("path",
             x = nose$x,
             y = nose$y,
             size = rel(1.5 * scaling)) +
    annotate("curve",
             x = ears$x[[1]],
             y = ears$y[[1]],
             xend = ears$x[[2]],
             yend = ears$y[[2]],
             curvature = -.5,
             angle = 60,
             size = rel(1.5 * scaling)) +
    annotate("curve",
             x = ears$x[[3]],
             y = ears$y[[3]],
             xend = ears$x[[4]],
             yend = ears$y[[4]],
             curvature = .5,
             angle = 120,
             size = rel(1.5 * scaling)) +
    coord_equal() +
    theme_bw() +
    theme(rect = element_blank(),
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()) +
    guides(fill = guide_colorbar(title = expression(paste("Amplitude (",
                                                          mu, "V)")),
                                 title.position = "right",
                                 barwidth = rel(1),
                                 barheight = rel(6),
                                 title.theme = element_text(angle = 270)))

  # Add electrode points or names -------------------
  if (chan_marker == "point") {
    topo <- topo +
      annotate("point",
               x = scaled_x,
               y = scaled_y,
               colour = "black",
               size = rel(2 * scaling))
    }  else if (chan_marker == "name") {
      topo <- topo +
        annotate("text",
                 x = scaled_x,
                 y = scaled_y,
                 label = c(levels(data$electrode)[c(data$electrode)]),
                 colour = "black",
                 size = rel(4 * scaling))
    }

  # Highlight specified electrodes
  if (!is.null(highlights)) {
    high_x <- scaled_x[data$electrode %in% highlights]
    high_y <- scaled_y[data$electrode %in% highlights]
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

#' @describeIn topoplot Topographical plotting of \code{eeg_data} objects.
#' @export

topoplot.eeg_data <- function(data, time_lim = NULL,
                              limits = NULL,
                              chanLocs = NULL,
                              method = "Biharmonic",
                              r = NULL,
                              grid_res = 67,
                              palette = "RdBu",
                              interp_limit = "skirt",
                              contour = TRUE,
                              chan_marker = "point",
                              quantity = "amplitude",
                              montage = NULL,
                              highlights = NULL,
                              scaling = 1,
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
           scaling = scaling)
}


#' @describeIn topoplot Topographical plotting of \code{eeg_epochs} objects.
#' @export

topoplot.eeg_epochs <- function(data,
                                time_lim = NULL,
                                limits = NULL,
                                chanLocs = NULL,
                                method = "Biharmonic",
                                r = NULL,
                                grid_res = 67,
                                palette = "RdBu",
                                interp_limit = "skirt",
                                contour = TRUE,
                                chan_marker = "point",
                                quantity = "amplitude",
                                montage = NULL,
                                highlights = NULL,
                                scaling = 1,
                                groups = NULL,
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
           groups = groups
           )
}


#' @param component Component to plot (numeric)
#' @describeIn topoplot Topographical plot for \code{eeg_ICA} objects
#' @export
topoplot.eeg_ICA <- function(data,
                             component,
                             time_lim = NULL,
                             limits = NULL,
                             chanLocs = NULL,
                             method = "Biharmonic",
                             r = NULL,
                             grid_res = 67,
                             palette = "RdBu",
                             interp_limit = "skirt",
                             contour = TRUE,
                             chan_marker = "point",
                             quantity = "amplitude",
                             montage = NULL,
                             colourmap,
                             highlights = NULL,
                             scaling = scaling,
                             ...) {
  if (missing(component)) {
    stop("Component number must be specified for eeg_ICA objects.")
  }

  chan_info <- data$chan_info
  data <- data.frame(amplitude = data$mixing_matrix[, component],
                     electrode = data$mixing_matrix$electrode)
  topoplot(data,
           chanLocs = chan_info,
           time_lim = NULL)

}

#' @param freq_range Range of frequencies to average over.
#' @describeIn topoplot Topographical plotting of \code{eeg_tfr} objects.
#' @export

topoplot.eeg_tfr <- function(data,
                             time_lim = NULL,
                             limits = NULL,
                             chanLocs = NULL,
                             method = "Biharmonic",
                             r = NULL,
                             grid_res = 67,
                             palette = "RdBu",
                             interp_limit = "skirt",
                             contour = TRUE,
                             chan_marker = "point",
                             quantity = "power",
                             montage = NULL,
                             highlights = NULL,
                             scaling = 1,
                             freq_range = NULL,
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
           passed = TRUE)
}

#' Create topoplot
#'
#' @noRd

make_topo <- function(data,
                      time_lim,
                      limits,
                      chanLocs,
                      method,
                      r,
                      grid_res,
                      palette,
                      interp_limit,
                      contour,
                      chan_marker,
                      quantity,
                      montage,
                      highlights) {


}

#' Set palette and limits for topoplot
#'
#' @param topo ggplot2 object produced by topoplot command
#' @param palette Requested palette
#' @param limits Limits of colour scale
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @keywords internal

set_palette <- function(topo, palette, limits = NULL) {

  if (palette %in% c("magma", "inferno", "plasma",
                  "viridis", "A", "B", "C", "D")) {

    topo <- topo +
      viridis::scale_fill_viridis(option = palette,
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

#' Fit a GAM smooth across scalp.
#'
#' @param data Data frame from which to generate predictions
#' @param scaled_x x coordinates rescaled
#' @param scaled_y y coordinates rescaled
#' @param grid_res resolution of output grid
#' @importFrom mgcv gam
#' @keywords internal

gam_topo <- function(data,
                     scaled_x,
                     scaled_y,
                     grid_res) {

  data$x <- scaled_x
  data$y <- scaled_y
  spline_smooth <- mgcv::gam(z ~ s(x,
                                   y,
                                   bs = "ts",
                                   k = 40),
                             data = data)
  out_df <- data.frame(expand.grid(x = seq(min(data$x) * 2,
                                           max(data$x) * 2,
                                           length = grid_res),
                                   y = seq(min(data$y) * 2,
                                           max(data$y) * 2,
                                           length = grid_res)))

  out_df$amplitude <-  stats::predict(spline_smooth,
                                      out_df,
                                      type = "response")
  out_df

}

#' Create a biharmonic smooth across the scalp
#'
#' @param data Data frame from which to generate predictions
#' @param scaled_x x coordinates rescaled
#' @param scaled_y y coordinates rescaled
#' @param grid_res resolution of output grid
#' @importFrom purrr map
#' @keywords internal

biharm_topo <- function(data,
                        scaled_x,
                        scaled_y,
                        grid_res) {

  # Create the interpolation grid --------------------------
  xo <- seq(-1.4, 1.4, length = grid_res)
  yo <- seq(-1.4, 1.4, length = grid_res)

  xo <- matrix(rep(xo, grid_res),
               nrow = grid_res,
               ncol = grid_res)
  yo <- t(matrix(rep(yo, grid_res),
                 nrow = grid_res,
                 ncol = grid_res))
  xy <- scaled_x + scaled_y * sqrt(as.complex(-1))
  d <- matrix(rep(xy, length(xy)),
              nrow = length(xy),
              ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d ^ 2) * (log(d) - 1) #Green's function
  diag(g) <- 0
  weights <- qr.solve(g, data$z)
  xy <- t(xy)

  outmat <-
    purrr::map(xo + sqrt(as.complex(-1)) * yo,
               function (x) (abs(x - xy) ^ 2) *
                 (log(abs(x - xy)) - 1) ) %>%
    rapply(function (x) ifelse(is.nan(x),
                               0,
                               x),
           how = "replace") %>%
    purrr::map_dbl(function (x) x %*% weights)

  dim(outmat) <- c(grid_res, grid_res)

  out_df <- data.frame(x = xo[, 1],
                       outmat)
  names(out_df)[1:length(yo[1, ]) + 1] <- yo[1, ]
  out_df <- tidyr::gather(out_df,
                          key = y,
                          value = amplitude,
                          -x,
                          convert = TRUE)
}

#' Create head shape for plotting
#'
#' Creates data frames for plotting a head with ears and nose.
#'
#' @param r Radius of head
#' @param circ_rads Circle outline in radians
#' @keywords internal

create_head <- function(r, circ_rads) {

  head_shape <- data.frame(x = r * cos(circ_rads),
                           y = r * sin(circ_rads))
  #define nose position relative to head_shape
  nose <- data.frame(x = c(head_shape$x[[23]],
                           head_shape$x[[26]],
                           head_shape$x[[29]]),
                     y = c(head_shape$y[[23]],
                           head_shape$y[[26]] * 1.1,
                           head_shape$y[[29]]))

  ears <- data.frame(x = c(head_shape$x[[4]],
                           head_shape$x[[97]],
                           head_shape$x[[48]],
                           head_shape$x[[55]]),
                     y = c(head_shape$y[[4]],
                           head_shape$y[[97]],
                           head_shape$y[[48]],
                           head_shape$y[[55]]))

  head_out <- list("head_shape" = head_shape,
                   "nose" = nose,
                   "ears" = ears)
  head_out
}

#' @noRd
round_topo <- function(.data, interp_limit, r, circ_rads) {

  if (identical(interp_limit, "skirt")) {
    .data$incircle <- sqrt(.data$x ^ 2 + .data$y ^ 2) < 1.125
    mask_ring <- data.frame(x = 1.126 * cos(circ_rads),
                            y = 1.126 * sin(circ_rads))
    } else {
      .data$incircle <- sqrt(.data$x ^ 2 + .data$y ^ 2) < (r * 1.02)
      mask_ring <- data.frame(x = r * 1.03 * cos(circ_rads),
                              y = r * 1.03 * sin(circ_rads))
    }
  topo_out <- list("out_df" = .data[.data$incircle, ],
                   "mask_ring" = mask_ring)
  topo_out
}


parse_for_topo <- function(.data,
                           time_lim) {

  if (!is.null(time_lim)) {
    .data <- select_times(.data, time_lim)
  }
  .data
}
