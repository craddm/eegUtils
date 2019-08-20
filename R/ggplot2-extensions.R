#' @importFrom ggplot2 fortify
#' @export
ggplot2::fortify

#' @importFrom tibble as_tibble
#' @export
fortify.eeg_epochs <- function(model,
                               data,
                               ...) {
  tibble::as_tibble(as.data.frame(model,
                                  long = TRUE,
                                  stringsAsFactors = FALSE))
}

#' @export
fortify.eeg_data <- function(model,
                               data,
                               ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}

#'@export
fortify.eeg_ICA <- function(model,
                            data,
                            ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}

#' @export
fortify.eeg_tfr <- function(model,
                            data,
                            ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}

#' @export
fortify.eeg_evoked <- function(model,
                               data,
                               ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}

#' StatScalpmap
#'
#' @export
#' @keywords internal

StatScalpmap <-
  ggplot2::ggproto("StatScalpmap",
                   Stat,
                   required_aes = c("x",
                                    "y",
                                    "fill"),
                   compute_group = function(data,
                                            scales,
                                            grid_res,
                                            interp_limit,
                                            method = "biharmonic") {

                     data <- aggregate(fill ~ x + y,
                                       data = data,
                                       FUN = mean)

                     if (identical(method, "biharmonic")) {
                       data <- biharmonic(data,
                                          grid_res = grid_res,
                                          interp_limit = interp_limit)
                     } else {
                       data <- fit_gam_topo(data,
                                        grid_res = grid_res,
                                        interp_limit = interp_limit)
                     }
                     data

                     }
)

#' Create an interpolated scalp surface
#'
#' \code{stat_scalpmap} creates an interpolated surface for an irregular set of
#' x-y coordinates, as is typically required for a topographical EEG plot. Since
#' the surface should be approximately round, the function attempts to blank out
#' portions of the surface that lay outside the area within the electrodes.
#'
#' @inheritParams ggplot2::geom_raster
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param grid_res Resolution of the interpolation grid. (Defaults to 100
#'   points).
#' @param interp_limit Topoplot with a "skirt" or inside the "head".
#' @param method "biharmonic" or "gam"
#' @family topoplot functions
#' @export
stat_scalpmap <- function(mapping = NULL,
                          data = NULL,
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE,
                          grid_res = 300,
                          interpolate = FALSE,
                          interp_limit = "skirt",
                          method = "biharmonic",
                          ...) {
  ggplot2::layer(
    stat = StatScalpmap,
    data = data,
    mapping = mapping,
    geom = GeomRaster,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  interpolate = interpolate,
                  grid_res = grid_res,
                  interp_limit = interp_limit,
                  method = method,
                  ...)
  )
}

StatConts <-
  ggplot2::ggproto("StatConts",
                   Stat,
                   required_aes = c("x",
                                    "y",
                                    "z"),
                   compute_group = function(data,
                                            scales,
                                            grid_res,
                                            interp_limit,
                                            method = "biharmonic") {

                     data <- aggregate(z ~ x + y,
                                       data = data,
                                       FUN = mean)
                     data <- rename(data, fill = "z")

                     if (identical(method, "biharmonic")) {
                       data <- biharmonic(data,
                                          grid_res = grid_res,
                                          interp_limit = interp_limit)
                     } else {
                       data <- fit_gam_topo(data,
                                            grid_res = grid_res,
                                            interp_limit = interp_limit)
                     }
                     data <- rename(data, z = "fill")
                     data
                   }
  )

#' Create a topographical plot
#'
#' \code{geom_topo()} creates a topographical plot as a \code{ggplot2} object.
#' This function automatically combines a number of distinct geom_* and stat_*
#' functions to create a default topographical scalp map. Since geom_raster does
#' not allow unevenly spaced grids, the function creates an interpolated surface.
#'
#' @examples
#' library(ggplot2)
#' ggplot(demo_epochs, aes(x = x, y = y, fill = amplitude)) + geom_topo()
#' @inheritParams ggplot2::geom_raster
#' @param interp_limit Topoplot with a "skirt" or inside the "head".
#' @param chan_markers Defaults to "point". Mark electrode positions with points
#'   or text.
#' @param chan_size Size for channel markers, if any.
#' @param head_size Size of the head shape.
#' @param grid_res Smoothness of the interpolation grid.
#' @param method "biharmonic" or ""gam".
#' @family topoplot functions
#' @export
geom_topo <- function(mapping = NULL,
                      data = NULL,
                      stat = "identity",
                      position = "identity",
                      show.legend = NA,
                      na.rm = TRUE,
                      inherit.aes = TRUE,
                      interpolate = FALSE,
                      interp_limit = "skirt",
                      chan_markers = "point",
                      chan_size = rel(2),
                      head_size = rel(1.5),
                      grid_res = 200,
                      method = "biharmonic",
                      ...) {

  list(ggplot2::layer(geom = GeomRaster,
                      stat = StatScalpmap,
                      data = data,
                      mapping = mapping,
                      position = position,
                      show.legend = show.legend,
                      inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm,
                                    interpolate = interpolate,
                                    grid_res = grid_res,
                                    interp_limit = interp_limit,
                                    method = method,
                                    ...)
                      ),
       ggplot2::layer(geom = GeomHead,
                      data = data,
                      mapping = mapping,
                      stat = StatHead,
                      position = PositionIdentity,
                      inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm,
                                    size = head_size,
                                    interp_limit = interp_limit,
                                    ...)
                      ),
       ggplot2::layer(data = data,
                      mapping = mapping,
                      stat = StatREar,
                      geom = GeomEars,
                      position = PositionIdentity,
                      show.legend = show.legend,
                      inherit.aes = TRUE,
                      params = list(na.rm = na.rm,
                                    curvature = -.5,
                                    angle = 60,
                                    size = head_size,
                                    interp_limit = interp_limit,
                                    ...)
                      ),
       ggplot2::layer(data = data,
                      mapping = mapping,
                      stat = StatLEar,
                      geom = GeomEars,
                      position = PositionIdentity,
                      show.legend = show.legend,
                      inherit.aes = TRUE,
                      params = list(na.rm = na.rm,
                                    curvature = .5,
                                    angle = 120,
                                    size = head_size,
                                    interp_limit = interp_limit,
                                    ...)
                      ),
       if (chan_markers == "point") {
         ggplot2::layer(data = data,
                        mapping = mapping,
                        stat = StatChannels,
                        geom = GeomPoint,
                        position = PositionIdentity,
                        show.legend = show.legend,
                        inherit.aes = inherit.aes,
                        params = list(na.rm = na.rm,
                                      fill = NA,
                                      size = chan_size,
                                      ...))
         } else if (chan_markers == "text") {
           ggplot2::layer(data = data,
                          mapping = mapping,
                          stat = StatChannels,
                          geom = GeomText,
                          position = PositionIdentity,
                          show.legend = show.legend,
                          inherit.aes = inherit.aes,
                          params = list(na.rm = na.rm,
                                        size = chan_size,
                                        ...))
           }
       )
}

#' @keywords internal
GeomTopo <- ggplot2::ggproto("GeomTopo",
                             GeomRaster)


#' Add head shape
#'
#' \code{geom_head()} adds a headshape to a plot.
#' @rdname stat_scalpmap
#' @inheritParams ggplot2::geom_path
#' @param r Radius of head
#' @family topoplot functions
#' @export
geom_head <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = TRUE,
                      inherit.aes = TRUE,
                      interp_limit = "skirt",
                      r = NULL,
                      ...) {

  list(ggplot2::layer(geom = GeomHead,
                      data = data,
                      mapping = mapping,
                      stat = StatHead,
                      position = PositionIdentity,
                      inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm,
                                    interp_limit = interp_limit,
                                    r = r,
                                    ...)),
       ggplot2::layer(data = data,
                      mapping = mapping,
                      stat = StatREar,
                      geom = GeomEars,
                      position = PositionIdentity,
                      show.legend = show.legend,
                      inherit.aes = TRUE,
                      params = list(na.rm = na.rm,
                                    curvature = -.5,
                                    angle = 60,
                                    interp_limit = interp_limit,
                                    r = r,
                                    ...)),
       ggplot2::layer(data = data,
                      mapping = mapping,
                      stat = StatLEar,
                      geom = GeomEars,
                      position = PositionIdentity,
                      show.legend = show.legend,
                      inherit.aes = TRUE,
                      params = list(na.rm = na.rm,
                                    curvature = .5,
                                    angle = 120,
                                    interp_limit = interp_limit,
                                    r = r,
                                    ...))
  )
}

StatHead <- ggplot2::ggproto("StatHead",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = NULL) {

                               if (!is.null(r)) {
                                 heads <- make_head(r = r)
                                 return(heads)
                               }

                               if (identical(interp_limit, "head")) {
                                 y_lim <- max(abs(data$y),
                                              na.rm = TRUE) * 1.1
                               } else {
                                 y_lim <- max(data$y,
                                              na.rm = TRUE) * 1.1
                               }

                               heads <- make_head(r = y_lim)
                               heads
                             }
)

GeomHead <- ggplot2::ggproto("GeomHead",
                             GeomPath)


#' Mask ring
#'
#' \code{geom_mask()} adds a masking ring to smooth the edges of a scalp map
#' generated by \code{stat_scalpmap()}, to give it a circular appearance.
#'
#' @inheritParams ggplot2::geom_path
#' @param colour For \code{geom_mask}, colour of the masking ring.
#' @param size For \code{geom_mask}, width of the masking ring.
#' @param scale_fac The radius of the ring is determined from the front-most
#'   electrode's location by a scaling factor. Defaults to 1.32 * max(y).
#' @rdname stat_scalpmap
#' @family topoplot functions
#' @export
geom_mask <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = FALSE,
                      colour = "white",
                      size = rel(6.5),
                      scale_fac = 1.32,
                      ...) {

  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatMask,
                 geom = GeomPath,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = TRUE,
                 params = list(na.rm = na.rm,
                               colour = colour,
                               size = size,
                               scale_fac = scale_fac,
                               ...))
}

StatMask <-
  ggplot2::ggproto("StatMask",
                   Stat,
                   compute_group = function(data,
                                            scales,
                                            scale_fac = 1.32) {

                     scale_fac <- scale_fac * max(data$y,
                                                  na.rm = TRUE)
                     data <- data.frame(x = scale_fac * cos(circ_rad_fun()),
                                        y = scale_fac * sin(circ_rad_fun()))
                     data

                   }
  )


#' Add ears to head
#'
#' \code{geom_ears} simply draws a pair of ears attached to the head shape.
#' @inheritParams ggplot2::geom_curve
#' @rdname stat_scalpmap
#' @export
geom_ears <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = FALSE,
                      ...) {

  list(
    ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatREar,
                 geom = GeomEars,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = TRUE,
                 params = list(na.rm = na.rm,
                               curvature = -.5,
                               angle = 60,
                               ...)),
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = StatLEar,
                   geom = GeomEars,
                   position = PositionIdentity,
                   show.legend = show.legend,
                   inherit.aes = TRUE,
                   params = list(na.rm = na.rm,
                                 curvature = .5,
                                 angle = 120,
                                 ...))
  )

}


GeomEars <- ggplot2::ggproto("GeomEars",
                             GeomCurve)

StatREar <- ggplot2::ggproto("StatREar",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = NULL) {

                               if (!is.null(r)) {
                                 return(make_r_ear(r = r))
                               }

                               if (identical(interp_limit, "head")) {
                                 y_lim <- max(abs(data$y),
                                              na.rm = TRUE) * 1.1
                               } else {
                                 y_lim <- max(data$y,
                                              na.rm = TRUE) * 1.1
                               }
                               make_r_ear(y_lim)
                             })

StatLEar <- ggplot2::ggproto("StatLEar",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = NULL) {

                               if (!is.null(r)) {
                                 return(make_l_ear(r = r))
                               }

                               if (identical(interp_limit, "head")) {
                                 y_lim <- max(abs(data$y),
                                              na.rm = TRUE) * 1.1
                               } else {
                                 y_lim <- max(data$y,
                                              na.rm = TRUE) * 1.1
                               }
                               make_l_ear(r = y_lim)
                             })

#' Create a headshape
#'
#' @keywords internal
make_head <- function(r) {

  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()),
                           group = 1)
  #define nose position relative to head_shape
  nose <- data.frame(x = c(head_shape$x[[23]],
                           head_shape$x[[26]],
                           head_shape$x[[29]]),
                     y = c(head_shape$y[[23]],
                           head_shape$y[[26]] * 1.1,
                           head_shape$y[[29]]),
                     group = 2)

  head_out <- rbind(head_shape, nose)
  head_out
}

#' Make right ear
#' @param r Radius of head
#' @keywords internal
make_r_ear <- function(r) {

  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()))
  right_ear <- data.frame(x = head_shape$x[[4]],
                          xend = head_shape$x[[97]],
                          y = head_shape$y[[4]],
                          yend = head_shape$y[[97]])
  right_ear
}

#' Make left ear
#' @param r Radius of head
#' @keywords internal
make_l_ear <- function(r) {
  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()))
  left_ear <- data.frame(x = head_shape$x[[48]],
                         xend = head_shape$x[[55]],
                         y = head_shape$y[[48]],
                         yend = head_shape$y[[55]])
  left_ear
}

StatChannels <-
  ggplot2::ggproto("StatChannels",
                   Stat,
                   required_aes = c("x", "y"),
                   compute_group = function(data, scales) {

                     if ("label" %in% names(data)) {
                     data <- aggregate(data[, c("x", "y")],
                                       by = list(label = data$label),
                                       FUN = mean)
                     } else {
                       data <- data[!duplicated(data[, c("x", "y")]),
                                    c("x", "y")]
                     }
                     })


#' Add channel indicators
#'
#' \code{geom_channels} adds either points or text labels at channel locations.
#' This is a convenience function to prevent overplotting when the input data
#' contains many rows of data.
#'
#' @inheritParams ggplot2::geom_point
#' @inheritParams ggplot2::geom_text
#' @rdname stat_scalpmap
#' @param geom "point" for points or "text" for labels. Default is "point".
#' @export
geom_channels <- function(mapping = NULL,
                          data = NULL,
                          geom = "point",
                          show.legend = NA,
                          inherit.aes = TRUE,
                          na.rm = TRUE,
                          ...) {

  ggplot2::layer(data = data,
                 mapping = mapping,
                 stat = StatChannels,
                 geom = geom,
                 position = PositionIdentity,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm,
                               ...))
}

biharmonic <- function(data,
                       grid_res,
                       interp_limit) {

  x_min <- min(data$x, na.rm = TRUE) * 2.5
  x_max <- max(data$x, na.rm = TRUE) * 2.5
  y_min <- min(data$y, na.rm = TRUE) * 2.5
  y_max <- max(data$y, na.rm = TRUE) * 2.5

  abs_y_max <- max(abs(data$y), na.rm = TRUE)
  xo <- seq(x_min,
            x_max,
            length = grid_res)
  yo <- seq(y_min,
            y_max,
            length = grid_res)
  xo <- matrix(rep(xo,
                   grid_res),
               nrow = grid_res,
               ncol = grid_res)
  yo <- t(matrix(rep(yo, grid_res),
                 nrow = grid_res,
                 ncol = grid_res))

  xy_coords <- unique(data[, c("x", "y")])
  xy <- xy_coords[, 1, drop = TRUE] + xy_coords[, 2, drop = TRUE] * sqrt(as.complex(-1))
  d <- matrix(rep(xy,
                  length(xy)),
              nrow = length(xy),
              ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d ^ 2) * (log(d) - 1) #Green's function
  diag(g) <- 0
  weights <- qr.solve(g, data$fill)
  xy <- t(xy)
  # Remind me to make this code readable at some point.
  outmat <-
    purrr::map(xo + sqrt(as.complex(-1)) * yo,
               function(x) (abs(x - xy) ^ 2) *
                 (log(abs(x - xy)) - 1)) %>%
    rapply(function(x) ifelse(is.nan(x), 0, x),
           how = "replace") %>%
    purrr::map_dbl(function(x) x %*% weights)

  dim(outmat) <- c(grid_res, grid_res)
  data <- data.frame(x = xo[, 1],
                     outmat)
  names(data)[1:length(yo[1, ]) + 1] <- yo[1, ]
  data <- tidyr::gather(data,
                        key = y,
                        value = fill,
                        -x,
                        convert = TRUE)

  if (identical(interp_limit, "head")) {
    circ_scale <- abs_y_max * 1.1
  } else {
    circ_scale <- y_max / 1.8
  }

  data[sqrt(data$x ^ 2 + data$y ^ 2) < circ_scale, ]
}

fit_gam_topo <- function(data,
                         grid_res,
                         interp_limit) {

  abs_y_max <- max(abs(data$y), na.rm = TRUE)
  y_max <- max(data$y, na.rm = TRUE) * 1.5

  spline_smooth <- mgcv::gam(fill ~ s(x,
                                      y,
                                      bs = "ts",
                                      k = 40),
                             data = data)
  data <- data.frame(expand.grid(x = seq(min(data$x) * 1.5,
                                         max(data$x) * 1.5,
                                         length = grid_res),
                                 y = seq(min(data$y) * 1.5,
                                         max(data$y) * 1.5,
                                         length = grid_res)))
  data$fill <-  stats::predict(spline_smooth,
                               data,
                               type = "response")


  if (identical(interp_limit, "head")) {
    circ_scale <- abs_y_max * 1.1
  } else {
    circ_scale <- y_max / 1.8
  }


  data$incircle <- sqrt(data$x ^ 2 + data$y ^ 2) < circ_scale
  data[data$incircle, ]
}


get_scalpmap <- function(data,
                         grid_res = 100,
                         interp_limit = "skirt",
                         method = "biharmonic") {

  # only useful if passed eeg objects
  # check_locs <- no_loc_chans(channels(data))
  # if (!is.null(check_locs)) {
  #   data <- select(data, -check_locs)
  # }
  # tmp <- as.data.frame(data)
  # tmp <- dplyr::summarise_at(tmp,
  #                            channel_names(data),
  #                            .funs = mean)
  # tmp <- tidyr::gather(tmp,
  #                      electrode,
  #                      amplitude)
  # tmp <- dplyr::left_join(tmp,
  #                         channels(data))
  tmp <- dplyr::filter(data, !is.na(x))
  tmp <- group_by(tmp,
                  electrode) %>%
    summarise(amplitude = mean(amplitude),
              x = mean(x), y= mean(y))
  smooth <-
    switch(method,
           biharmonic = biharmonic(dplyr::rename(tmp,
                                                 fill = "amplitude"),
                                   grid_res = grid_res,
                                   interp_limit = interp_limit),
           gam = fit_gam_topo(dplyr::rename(tmp,
                                            fill = "amplitude"),
                              grid_res = grid_res,
                              interp_limit = interp_limit))

  smooth
}


get_scalpmapn <- function(data, ...) {
  UseMethod("getscalpmapn", data)
}

get_scalpmapn.data.frame <- function(data, ...) {

}


get_scalpmapn.eeg_epochs <- function(data, ...) {
  # only useful if passed eeg objects
   check_locs <- no_loc_chans(channels(data))
   if (!is.null(check_locs)) {
     data <- select(data, -check_locs)
   }
  tmp <- as.data.frame(data)
  tmp <- dplyr::summarise_at(tmp,
                             channel_names(data),
                             .funs = mean)
  tmp <- tidyr::gather(tmp,
                       electrode,
                       amplitude)
  tmp <- dplyr::left_join(tmp,
                          channels(data))
}
no_loc_chans <- function(chaninfo) {

  if (any(is.na(chaninfo$x))) {
    return(chaninfo$electrode[is.na(chaninfo$x)])
  }
  NULL
}
