#' Topographical Plotting Function for EEG
#'
#' Allows simple plotting of functional data. Output is a ggplot2 object.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @param data An EEG dataset. Must have columns x, y, and amplitude at present.
#'   x and y are (Cartesian) electrode co-ordinates), amplitude is amplitude.
#' @param time_lim Timepoint(s) to plot. Can be one time or a range to average
#'   over. If none is supplied, the function will average across all timepoints
#'   in the supplied data.
#' @param clim Limits of the fill scale - should be given as a character vector
#'   with two values specifying the start and endpoints e.g. clim = c(-2,-2).
#'   Will ignore anything else. Defaults to the range of the data.
#' @param chanLocs Not yet implemented.
#' @param method Interpolation method. "Biharmonic" or "gam". "Biharmonic"
#'   implements the same method used in Matlab's EEGLAB. "gam" fits a
#'   Generalized Additive Model with k = 40 knots. Defaults to biharmonic spline
#'   interpolation.
#' @param r Radius of cartoon headshape; if not given, defaults to 1.1 * the
#'   maximum y electrode location.
#' @param gridRes Resolution of the interpolated grid. Higher = smoother but
#'   slower.
#' @param colourmap Defaults to RdBu if none supplied. Can be any from
#'   RColorBrewer. If an unsupported palette is specified, switches to Greens.
#' @param interp_limit "skirt" or "head". Defaults to "skirt". "skirt"
#'   interpolates just past the farthest electrode and does not respect the
#'   boundary of the headshape. "head" interpolates up to the radius of the
#'   plotted head.
#' @param contour Plot contour lines on topography (defaults to TRUE)
#' @param chan_marker Set marker for electrode locations. "point" = point,
#'   "name" = electrode name, "none" = no marker. Defaults to "point".
#' @param quantity Allows plotting of arbitrary quantitative column. Defaults to
#'   amplitude. Can be any column name. E.g. "p.value", "t-statistic".
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom rlang parse_quosure
#' @import scales
#' @importFrom mgcv gam
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

topoplot <- function(data,
                     time_lim = NULL,
                     clim = NULL,
                     chanLocs = NULL,
                     method = "Biharmonic",
                     r = NULL,
                     gridRes = 67,
                     colourmap = "RdBu",
                     interp_limit = "skirt",
                     contour = TRUE,
                     chan_marker = "point",
                     quantity = "amplitude",
                     montage = NULL) {

  # Check if data is of class eeg_data.
  if (is.eeg_data(data)) {
    data <- as.data.frame(data, long = TRUE)
  }

  # Filter out unwanted timepoints, and find nearest time values in the data --------------

  if ("time" %in% colnames(data)) {
    if (length(time_lim) == 1) {
      time_lim <- data$time[which.min(abs(data$time - time_lim))]
      data <- dplyr::filter(data, time == time_lim)
      } else if (length(time_lim) == 2) {
        data <- select_times(data, time_lim)
      }
    }

  # Check for x and y co-ordinates, try to add if not found --------------

  if (length(grep("^x$|^y$", colnames(data))) > 1) {
    message("Electrode locations found.")
  } else if (!is.null(chanLocs)) {
    if (length(grep("^x$|^y$", colnames(chanLocs))) > 1) {
      data <- dplyr::left_join(data, chanLocs, by = "electrode")
      } else {
        warnings("No channel locations found in chanLocs.")
      }
    } else if ("electrode" %in% colnames(data)) {
      data <- electrode_locations(data, drop = TRUE, montage = montage)
      message("Attempting to add standard electrode locations...")
    } else {
    warning("Neither electrode locations nor labels found.")
    stop()
  }

  # Average over all timepoints ----------------------------

  data <- dplyr::summarise(dplyr::group_by(data, x, y, electrode),
                   z = mean(!!rlang::parse_quosure(quantity)))

  # Cut the data frame down to only the necessary columns, and make sure it has the right names --------
  # Will be able to work with different parameters (e.g. power) eventually, so this will be necessary
  data <- data.frame(x = data$x,
                   y = data$y,
                   z = data$z,
                   electrode = data$electrode)

  # Rescale electrode co-ordinates to be from -1 to 1 for plotting
  # Selects largest absolute value from x or y
  max_dim <- max(abs(data$x), abs(data$y))
  scaled_x <- data$x/max_dim
  scaled_y <- data$y/max_dim

  # Create the interpolation grid --------------------------

  xo <- seq(-1.4, 1.4, length = gridRes)
  yo <- seq(-1.4, 1.4, length = gridRes)

  # Create the headshape -----------------

  #set radius as max of y (i.e. furthest forward electrode's y position). Add a little to push the circle out a bit more.

  if (is.null(r)) {
    r <- max(scaled_y) * 1.1
  }

  circ_rads <- seq(0, 2 * pi, length.out = 101)

  headShape <- data.frame(x = r * cos(circ_rads),
                          y = r * sin(circ_rads))

  #define nose position relative to headShape
  nose <- data.frame(x = c(headShape$x[[23]],
                           headShape$x[[26]],
                           headShape$x[[29]]),
                     y = c(headShape$y[[23]],
                           headShape$y[[26]] * 1.1,
                           headShape$y[[29]]))

  ears <- data.frame(x = c(headShape$x[[4]],
                           headShape$x[[97]],
                           headShape$x[[48]],
                           headShape$x[[55]]),
                     y = c(headShape$y[[4]],
                           headShape$y[[97]],
                           headShape$y[[48]],
                           headShape$y[[55]]))

  # Do the interpolation! ------------------------

  switch(method,
                  Biharmonic = {
                    xo <- matrix(rep(xo, gridRes), nrow = gridRes, ncol = gridRes)
                    yo <- t(matrix(rep(yo, gridRes), nrow = gridRes, ncol = gridRes))
                    xy <- scaled_x + scaled_y * sqrt(as.complex(-1))
                    d <- matrix(rep(xy, length(xy)), nrow = length(xy), ncol = length(xy))
                    d <- abs(d - t(d))
                    diag(d) <- 1
                    g <- (d ^ 2) * (log(d) - 1) #Green's function
                    diag(g) <- 0
                    weights <- qr.solve(g, data$z)
                    xy <- t(xy)

                    outmat <- purrr::map(xo + sqrt(as.complex(-1)) * yo,
                                      function (x) (abs(x - xy) ^ 2) * (log(abs(x - xy))-1) ) %>%
                      rapply(function (x) ifelse(is.nan(x), 0, x), how = "replace") %>%
                      purrr::map_dbl(function (x) x %*% weights)

                    dim(outmat) <- c(gridRes, gridRes)

                    outDf <- data.frame(x = xo[, 1], outmat)
                    names(outDf)[1:length(yo[1, ]) + 1] <- yo[1, ]
                    outDf <- tidyr::gather(outDf,
                                    key = y,
                                    value = amplitude,
                                    -x,
                                    convert = TRUE)
                  },
         gam = {
           tmp_df <- data
           tmp_df$x <- scaled_x
           tmp_df$y <- scaled_y
           splineSmooth <- mgcv::gam(z ~ s(x, y, bs = 'ts', k = 40), data = tmp_df)
           outDf <- data.frame(expand.grid(x = seq(min(tmp_df$x) * 2,
                                                     max(tmp_df$x) * 2,
                                                     length = gridRes),
                                             y = seq(min(tmp_df$y) * 2,
                                                     max(tmp_df$y) * 2,
                                                     length = gridRes)))

           outDf$amplitude <-  stats::predict(splineSmooth,
                                         outDf,
                                         type = "response")
         })

  # Check if should interp/extrap beyond headshape, and set up ring to mask edges for smoothness
  if (interp_limit == "skirt") {
    outDf$incircle <- sqrt(outDf$x ^ 2 + outDf$y ^ 2) < 1.125
    maskRing <- data.frame(x = 1.126 * cos(circ_rads),
                           y = 1.126 * sin(circ_rads)
    )
  } else {
    outDf$incircle <- sqrt(outDf$x ^ 2 + outDf$y ^ 2) < (r*1.03)
    maskRing <- data.frame(x = r * 1.03 * cos(circ_rads),
                           y = r * 1.03 * sin(circ_rads)
    )
  }

  # Create the actual plot -------------------------------

  topo <- ggplot2::ggplot(outDf[outDf$incircle, ], aes(x, y, fill = amplitude)) +
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
             x = maskRing$x,
             y = maskRing$y,
             colour = "white",
             size = rel(4.4)) +
    annotate("path",
             x = headShape$x,
             y = headShape$y,
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
    guides(fill = guide_colorbar(title = expression(paste("Amplitude (", mu,"V)")),
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

  # Set the colourmap and scale limits ------------------------

  if (length(clim) == 2) {
    topo <- topo + scale_fill_distiller(palette = colourmap,
                                limits = clim,
                                guide = "colourbar",
                                oob = scales::squish)
  } else {
    topo <- topo + scale_fill_distiller(palette = colourmap,
                                guide = "colourbar",
                                oob = scales::squish)
  }
  topo
}
