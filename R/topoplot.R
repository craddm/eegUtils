#' Topographical Plotting Function for EEG
#'
#' Allows topographical plotting of functional data. Output is a ggplot2 object.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data An EEG dataset. Must have columns x, y, and amplitude at present.
#'   x and y are (Cartesian) electrode co-ordinates), amplitude is amplitude.
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

#' Topographical Plotting Function for EEG
#'
#' @param data An object passed to the function
#' @param ... Any other parameters
#' @export
#'
topoplot.default <- function(data, ...) {
  stop("This function requires a data frame or an eeg_data/eeg_epochs object")
}

#' Topographical Plotting Function for EEG
#'
#' The functions works for both standard data frames and objects of class
#' \code{eeg_data}.
#'
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
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom rlang parse_quosure
#' @import scales
#' @importFrom mgcv gam
#' @describeIn topoplot Topographical plotting of data.frames and other non
#'   eeg_data objects.
#' @export


topoplot.data.frame <- function(data, time_lim = NULL, limits = NULL,
                             chanLocs = NULL, method = "Biharmonic", r = NULL,
                             grid_res = 67, palette = "RdBu",
                             interp_limit = "skirt", contour = TRUE,
                             chan_marker = "point", quantity = "amplitude",
                             montage = NULL, colourmap, ...) {

  if (!missing(colourmap)) {
    warning("Argument colourmap is deprecated, please use palette instead.",
            call. = FALSE)
    palette <- colourmap
  }

  # Filter out unwanted timepoints, and find nearest time values in the data
  # --------------

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      time_lim <- data$time[which.min(abs(data$time - time_lim))]
      data <- data[data$time == time_lim, ]
      } else if (length(time_lim) == 2) {
        range(data$time)
        data <- select_times(data, time_lim)
      }
    }

  # Check for x and y co-ordinates, try to add if not found --------------

  if (length(grep("^x$|^y$", colnames(data))) > 1) {
    message("Electrode locations found.")
  } else if (!is.null(chanLocs)) {
    if (length(grep("^x$|^y$", colnames(chanLocs))) > 1) {
      data$electrode <- toupper(data$electrode)
      data <- dplyr::left_join(data, chanLocs, by = "electrode")
      if (any(is.na(data$x))) {
        data <- data[!is.na(data$x),]
        }
      } else {
        warnings("No channel locations found in chanLocs.")
      }
    } else if ("electrode" %in% colnames(data)) {
      data <- electrode_locations(data, drop = TRUE, montage = montage)
      message("Attempting to add standard electrode locations...")
    } else {
    stop("Neither electrode locations nor labels found.")
  }

  # Average over all timepoints ----------------------------

  data <- dplyr::summarise(dplyr::group_by(data, x, y, electrode),
                           z = mean(!!rlang::parse_quosure(quantity)))

  # Cut the data frame down to only the necessary columns, and make sure it has
  # the right names
  data <- data.frame(x = data$x,
                   y = data$y,
                   z = data$z,
                   electrode = data$electrode)

  # Rescale electrode co-ordinates to be from -1 to 1 for plotting
  # Selects largest absolute value from x or y
  max_dim <- max(abs(data$x), abs(data$y))
  scaled_x <- data$x / max_dim
  scaled_y <- data$y / max_dim

  # Create the interpolation grid --------------------------

  xo <- seq(-1.4, 1.4, length = grid_res)
  yo <- seq(-1.4, 1.4, length = grid_res)

  # Create the head_shape -----------------

  #set radius as max of y (i.e. furthest forward electrode's y position). Add a
  #little to push the circle out a bit more.

  if (is.null(r)) {
    r <- max(scaled_y) * 1.1
  }

  circ_rads <- seq(0, 2 * pi, length.out = 101)

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

  # Do the interpolation! ------------------------

  switch(method,
         Biharmonic = {
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
             rapply(function (x) ifelse(is.nan(x), 0, x), how = "replace") %>%
             purrr::map_dbl(function (x) x %*% weights)

           dim(outmat) <- c(grid_res, grid_res)

           out_df <- data.frame(x = xo[, 1], outmat)
           names(out_df)[1:length(yo[1, ]) + 1] <- yo[1, ]
           out_df <- tidyr::gather(out_df,
                                  key = y,
                                  value = amplitude,
                                  -x,
                                  convert = TRUE)
         },
         gam = {
           out_df <- gam_topo(data, scaled_x, scaled_y, grid_res, 40)
         })

  # Check if should interp/extrap beyond head_shape, and set up ring to mask
  # edges for smoothness
  if (identical(interp_limit, "skirt")) {
    out_df$incircle <- sqrt(out_df$x ^ 2 + out_df$y ^ 2) < 1.125
    mask_ring <- data.frame(x = 1.126 * cos(circ_rads),
                           y = 1.126 * sin(circ_rads)
    )
  } else {
    out_df$incircle <- sqrt(out_df$x ^ 2 + out_df$y ^ 2) < (r * 1.03)
    mask_ring <- data.frame(x = r * 1.03 * cos(circ_rads),
                           y = r * 1.03 * sin(circ_rads)
    )
  }

  # Create the actual plot -------------------------------

  topo <- ggplot2::ggplot(out_df[out_df$incircle, ],
                          aes(x, y, fill = amplitude)) +
    geom_raster(interpolate = TRUE)

  if (contour) {
    topo <- topo + stat_contour(
      aes(z = amplitude, linetype = ..level.. < 0),
      bins = 6,
      colour = "black",
      size = 1.1,
      show.legend = FALSE
    )
  }

  topo <- topo +
    annotate("path",
             x = mask_ring$x,
             y = mask_ring$y,
             colour = "white",
             size = rel(4.4)) +
    annotate("path",
             x = head_shape$x,
             y = head_shape$y,
              size = rel(1.5)) +
    annotate("path",
             x = nose$x,
             y = nose$y,
              size = rel(1.5)) +
    annotate("curve",
             x = ears$x[[1]],
             y = ears$y[[1]],
             xend = ears$x[[2]],
             yend = ears$y[[2]],
             curvature = -.5,
             angle = 60,
             size = rel(1.5)) +
    annotate("curve",
             x = ears$x[[3]],
             y = ears$y[[3]],
             xend = ears$x[[4]],
             yend = ears$y[[4]],
             curvature = .5,
             angle = 120,
             size = rel(1.5)) +
    coord_equal() +
    theme_bw() +
    theme(rect = element_blank(),
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()) +
    guides(fill = guide_colorbar(title = expression(paste("Amplitude (",
                                                          mu, "V)")),
                                 title.position = "right",
                                 barwidth = 1,
                                 barheight = 6,
                                 title.theme = element_text(angle = 270)))

  # Add electrode points or names -------------------
  if (chan_marker == "point") {
    topo <- topo +
      annotate("point",
               x = scaled_x, y = scaled_y,
               colour = "black",
               size = rel(2))
    }  else if (chan_marker == "name") {
      topo <- topo +
        annotate("text",
                 x = scaled_x, y = scaled_y,
                 label = c(levels(data$electrode)[c(data$electrode)]),
                 colour = "black",
                 size = rel(4))
    }

  # Set the palette and scale limits ------------------------
  topo <- set_palette(topo, palette, limits)
  topo
}


#' Topographical Plotting Function for EEG data
#'
#' Both \code{eeg_epochs} and \code{eeg_data} objects are supported.
#'
#' @describeIn topoplot Topographical plotting of \code{eeg_data} objects.
#' @export

topoplot.eeg_data <- function(data, time_lim = NULL, limits = NULL,
                              chanLocs = NULL, method = "Biharmonic", r = NULL,
                              grid_res = 67, palette = "RdBu",
                              interp_limit = "skirt", contour = TRUE,
                              chan_marker = "point", quantity = "amplitude",
                              montage = NULL, ...) {

  if (!is.null(data$chan_info)) {
    chanLocs <- data$chan_info
  }
  data <- as.data.frame(data, long = TRUE)
  topoplot(data, time_lim = time_lim, limits = limits,
           chanLocs = chanLocs, method = method, r = r,
           grid_res = grid_res, palette = palette,
           interp_limit = interp_limit, contour = contour,
           chan_marker = chan_marker, quantity = quantity,
           montage = montage, passed = TRUE)
}


#' @param comp Component to plot (numeric)
#' @describeIn topoplot Topographical plot for \code{eeg_ICA} objects
topoplot.eeg_ICA <- function(data, time_lim = NULL, limits = NULL,
                             chanLocs = NULL, method = "Biharmonic", r = NULL,
                             grid_res = 67, palette = "RdBu",
                             interp_limit = "skirt", contour = TRUE,
                             chan_marker = "point", quantity = "amplitude",
                             montage = NULL, colourmap, comp, ...) {
  data <- data.frame(amplitude = scale(data$mixing_matrix[, comp]),
                     electrode = data$mixing_matrix$electrode)
  topoplot(data)

}

#' Set palette and limits for topoplot
#'
#' @param topo ggplot2 object produced by topoplot command
#' @param palette Requested palette
#' @param limits Limits of colour scale
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis


set_palette <- function(topo, palette, limits = NULL) {

  if (palette %in% c("magma", "inferno", "plasma",
                  "viridis", "A", "B", "C", "D")) {

    topo <- topo + viridis::scale_fill_viridis(option = palette,
                                      limits = limits,
                                      guide = "colourbar",
                                      oob = scales::squish)
  } else {
    topo <- topo + scale_fill_distiller(palette = palette,
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
#' @param k number of knots for GAM

gam_topo <- function(data, scaled_x, scaled_y, grid_res, k = 40) {

  data$x <- scaled_x
  data$y <- scaled_y
  spline_smooth <- mgcv::gam(z ~ s(x, y, bs = "ts", k = 40),
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
