#' Topographical Plotting Function for EEG
#'
#' Allows simple plotting of functional data. Output is a ggplot2 object.
#'
#' @param df An EEG dataset. Must have columns x, y, and amplitude at present. x and y are (Cartesian) electrode co-ordinates), amplitude is amplitude.
#' @param timepoint Timepoint(s) to plot. Can be one timepoint or a range to average over. If none is supplied, the function will average across all timepoints in the supplied data.
#' @param clims Limits of the fill scale - should be given as a character vector with two values specifying the start and endpoints e.g. clim = c(-2,-2). Will ignore anything else. Defaults to the range of the data.
#' @param chanLocs Not yet implemented.
#' @param method Interpolation method. "Biharmonic" or "gam". "Biharmonic" implements the same method used in Matlab's EEGLAB. "gam" fits a Generalized Additive Model Defaults to biharmonic spline interpolation.
#' @param rmax Maximum head radius if x is too small.
#' @param gridRes Resolution of the interpolated grid. Higher = smoother but slower.
#' @param colourmap Defaults to RdBu if none supplied. Can be any from RColorBrewer. If an unsupported palette is specified, switches to Greens.
#' @param chanNames Switch from points to electrode names if available. Not yet implemented.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @section Notes on usage of Generalized Additive Models for interpolation:
#' The function fits a GAM using the gam function from mgcv. Specifically, it fits a spline using the model function gam(z, s(x, y, bs = "sos"). Using GAMs for smooths is very much experimental. The surface is produced from the predictions of the GAM model fitted to the supplied data. Values at each electrode do not necessarily match actual values in the data: high-frequency variation will tend to be smoothed out. Thus, the method should be used with caution.

topoplot <- function(df,
                     timepoint = NULL,
                     clim = NULL,
                     chanLocs = NULL,
                     method = "Biharmonic",
                     rmax = .75,
                     gridRes = 67,
                     colourmap = "RdBu",
                     chanNames = FALSE) {

  # Filter out unwanted timepoints, and find nearest time values in the data --------------

  if ("time" %in% colnames(df)) {
    if (length(timepoint) == 1) {
      timepoint <- df$time[which.min(abs(df$time - timepoint))]
      df <- filter(df, time == timepoint)
      } else if (length(timepoint) == 2) {
        timepoint[1] <- df$time[which.min(abs(df$time - timepoint[1]))]
        timepoint[2] <- df$time[which.min(abs(df$time - timepoint[2]))]
        df <- filter(df, time >= timepoint[1] & time <= timepoint[2])
      }
    }

  # Check for x and y co-ordinates, try to add if not found --------------

  if (length(grep("^x$|^y$", colnames(df))) > 1) {
    message("Electrode locations found.")
  } else if ("electrode" %in% colnames(df)) {
    df <- left_join(df, electrodeLocs, by = "electrode")
    message("Adding standard electrode locations...")
  } else {
    warning("Neither electrode locations nor labels found.")
    stop()
  }

  # Average over all timepoints ----------------------------

  df <- summarise(group_by(df, x, y, electrode), amplitude = mean(amplitude))

  # Cut the data frame down to only the necessary columns, and make sure it has the right names --------
  # Will be able to work with different parameters (e.g. power) eventually, so this will be necessary
  df <- data.frame(x = df$x, y = df$y, z = df$amplitude, electrode = df$electrode)

  # Create the interpolation grid --------------------------

  xo <- seq(min(-rmax, df$x), max(rmax, df$x), length = gridRes)
  yo <- seq(max(rmax, df$y), min(-rmax, df$y), length = gridRes)

  # Create the headshape -----------------

  r <- round(max(df$x)) / 2;
  tt <- seq(0, 2 * pi, length.out = 100)
  headShape <- data.frame(x = r * cos(tt),
                           y = r * sin(tt))

  maskRing <- data.frame(x = (1.42/2) * cos(tt),
                         y = (1.42/2) * sin(tt)
                         )

  nose <- data.frame(x = c(-0.05, 0, .05), y = c(.495, .55, .495))

  # Do the interpolation! ------------------------

  switch(method,
                  Biharmonic = {
                    xo <- matrix(rep(xo, gridRes), nrow = gridRes, ncol = gridRes)
                    yo <- t(matrix(rep(yo, gridRes), nrow = gridRes, ncol = gridRes))
                    xy <- df$x + df$y * sqrt(as.complex(-1))
                    d <- matrix(rep(xy, length(xy)), nrow = length(xy), ncol = length(xy))
                    d <- abs(d - t(d))
                    diag(d) <- 1
                    g <- (d ^ 2) * (log(d) - 1) #Green's function
                    diag(g) <- 0
                    weights <- qr.solve(g, df$z)
                    xy <- t(xy)

                    outmat <- purrr::map(xo + sqrt(as.complex(-1)) * yo,
                                      function (x) (abs(x - xy) ^ 2) * (log(abs(x - xy))-1) ) %>%
                      rapply(function (x) ifelse(is.nan(x),0,x), how = "replace") %>%
                      purrr::map_dbl(function (x) x %*% weights)

                    dim(outmat) <- c(67,67)

                    outDf <- data.frame(x = xo[, 1], outmat)
                    names(outDf)[1:length(yo[1, ]) + 1] <- yo[1, ]
                    outDf <- gather(outDf,
                                    key = y,
                                    value = amplitude,
                                    -x,
                                    convert = TRUE)
                  },
         gam = {
           splineSmooth <- mgcv::gam(z ~ s(x, y, bs = 'sos'), data = df)
           outDf <- data.frame(expand.grid(x = seq(min(df$x) * 2,
                                                     max(df$x) * 2,
                                                     length = gridRes),
                                             y = seq(min(df$y) * 2,
                                                     max(df$y) * 2,
                                                     length = gridRes)))

           outDf$amplitude <-  predict(splineSmooth,
                                         outDf,
                                         type = "response")
         })

  outDf$incircle <- sqrt(outDf$x ^ 2 + outDf$y ^ 2) < .72

  # Create the actual plot -------------------------------

  topo <- ggplot(outDf[outDf$incircle, ], aes(x, y, fill = amplitude)) +
    geom_raster(interpolate = TRUE)+
    stat_contour(aes(z = amplitude, linetype = ..level..<0),
                 bins = 6,
                 colour = "black",
                 size = 1.2,
                 show.legend = FALSE) +
    annotate("path",
             x = headShape$x,
             y = headShape$y,
              size = 1.5) +
    annotate("path",
             x = nose$x,
             y = nose$y,
              size = 1.5) +
    annotate("path",
             x = maskRing$x,
             y = maskRing$y,
             colour = "white",
             size = 6)+
    coord_equal() +
    theme_bw() +
    theme(rect = element_blank(),
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()) +
    guides(fill = guide_colorbar(title = "Amplitude (microV)",
                                  title.position = "right",
                                  barwidth = 1,
                                  barheight = 6,
                                  title.theme = element_text(angle = 270)))

  # Add electrode points or names -------------------
  if (chanNames == TRUE) {
    topo <- topo + annotate("text",
                            x = df$x,
                            y = df$y,
                            label = c(levels(df$electrode)[c(df$electrode)]),
                            colour = "black",
                            size = 4.5)
  }
  else {
    topo <- topo + annotate("point",
           x = df$x,
           y = df$y,
           colour = "black",
           size = 2)
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
