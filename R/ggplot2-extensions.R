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

#' Create an interpolated scalp surface
#'
#' `stat_scalpmap` creates an interpolated surface for an irregular set of
#' x-y coordinates, as is typically required for a topographical EEG plot. Since
#' the surface should be approximately round, the function attempts to blank out
#' portions of the surface that lay outside the area within the electrodes.
#'
#' @inheritParams ggplot2::geom_raster
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param grid_res Resolution of the interpolation grid. (Defaults to 200
#'   points).
#' @param interp_limit Topoplot with a "skirt" or inside the "head".
#' @param method "biharmonic" or "gam"
#' @param r Size of head
#' @family topoplot functions
#' @export
stat_scalpmap <- function(mapping = NULL,
                          data = NULL,
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE,
                          grid_res = 200,
                          interpolate = FALSE,
                          interp_limit = c("skirt", "head"),
                          method = "biharmonic",
                          r = NULL,
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
                  r = r,
                  ...)
  )
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
                                            method,
                                            r) {

                     data <- aggregate(fill ~ x + y,
                                       data = data,
                                       FUN = mean)

                     # Add head and mask to topoplot

                     if (is.null(r)) {
                       abs_x_max <- max(abs(data$x), na.rm = TRUE)
                       abs_y_max <- max(abs(data$y), na.rm = TRUE)
                       r <- switch(interp_limit,
                                   "head" = sqrt(abs_x_max^2 + abs_y_max^2),
                                   "skirt" = 95) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
                     }

                     if (identical(method, "biharmonic")) {
                       data <- biharmonic(data,
                                          grid_res = grid_res,
                                          interp_limit = interp_limit,
                                          r = r)
                     } else {
                       data <- fit_gam_topo(data,
                                            grid_res = grid_res,
                                            interp_limit = interp_limit,
                                            r = r)
                     }
                     data

                   }
  )

#' Create a topographical plot
#'
#' `geom_topo()` creates a topographical plot as a `ggplot2` object.
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
#' @param r Head circumference
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
                      r = NULL,
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
                                    r = r,
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
                                    r = r,
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
                                    r = r,
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
                                    r = r,
                                    interp_limit = interp_limit,
                                    ...)
                      ),
       if (identical(chan_markers,
                     "point")) {
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
         } else if (identical(chan_markers,
                              "text")) {
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
#' `geom_head()` adds a headshape to a plot.
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
                      r = 95,
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
                                                      r = 95) {
                               # Add head and mask to topoplot
                               if (is.null(r)) {
                                 abs_x_max <- max(abs(data$x),
                                                  na.rm = TRUE)
                                 abs_y_max <- max(abs(data$y),
                                                  na.rm = TRUE)
                                 r <- switch(interp_limit,
                                             "head" = sqrt(abs_x_max^2 + abs_y_max^2),
                                             "skirt" = 95) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
                               }
                               r <- update_r(r,
                                             data,
                                             interp_limit)
                               make_head(r = r)
                             }
)

GeomHead <- ggplot2::ggproto("GeomHead",
                             GeomPath)


#' Mask ring
#'
#' `geom_mask()` adds a masking ring to smooth the edges of a scalp map
#' generated by `stat_scalpmap()`, to give it a circular appearance.
#'
#' @inheritParams ggplot2::geom_path
#' @param colour For `geom_mask`, colour of the masking ring.
#' @param size For `geom_mask`, width of the masking ring.
#' @param scale_fac The radius of the ring is determined from the front-most
#'   electrode's location by a scaling factor. Defaults to 1.04 * max(y).
#' @rdname stat_scalpmap
#' @family topoplot functions
#' @export
geom_mask <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = FALSE,
                      colour = "white",
                      size = rel(5),
                      scale_fac = 1.04,
                      r = 95,
                      interp_limit = "skirt",
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
                               r = r,
                               interp_limit = interp_limit,
                               ...))
}

StatMask <-
  ggplot2::ggproto("StatMask",
                   Stat,
                   compute_group = function(data,
                                            scales,
                                            scale_fac,
                                            interp_limit,
                                            r) {


                     abs_x_max <- max(abs(data$x),
                                      na.rm = TRUE)
                     abs_y_max <- max(abs(data$y),
                                      na.rm = TRUE)

                     scale_fac <- max(abs_x_max,
                                      abs_y_max) * scale_fac

                     data <- data.frame(x = scale_fac * cos(circ_rad_fun()),
                                        y = scale_fac * sin(circ_rad_fun()))
                     data

                   }
  )


#' Add ears to head
#'
#' `geom_ears` simply draws a pair of ears attached to the head shape.
#' @inheritParams ggplot2::geom_curve
#' @rdname stat_scalpmap
#' @export
geom_ears <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = FALSE,
                      r = 95,
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
                                 r = r,
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
                                                      r = 95) {

                               if (is.null(r)) {
                                 r <- 95
                               }
                               r <- update_r(r,
                                             data,
                                             interp_limit)
                               make_r_ear(r = r)
                             })

StatLEar <- ggplot2::ggproto("StatLEar",
                             Stat,
                             compute_group = function(data,
                                                      scales,
                                                      interp_limit,
                                                      r = NULL) {
                               if (is.null(r)) {
                                 r <- 95
                               }
                               r <- update_r(r,
                                             data,
                                             interp_limit)
                               make_l_ear(r = r)
                             }
)

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

  head_out <- rbind(head_shape,
                    nose)
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
#' `geom_channels` adds either points or text labels at channel locations.
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

stat_summary_by_fill <- function(mapping = NULL,
                                 data = NULL,
                                 geom = "raster",
                                 position = "identity",
                                 fun.data = NULL,
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = TRUE,
                                 ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatSummaryByFill,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      fun.data = fun.data,
      na.rm = na.rm,
      ...
    )
  )
}


StatSummaryByFill <- ggproto("StatSummaryByFill",
                             Stat,
                             required_aes = c("x", "y", "fill"),
                             compute_group = function(data,
                                                      scales,
                                                      fun.data = NULL,
                                                      na.rm = FALSE,
                                                      params,
                                                      layout) {
                                summary <-
                                  aggregate(fill ~ x + y,
                                            data = data,
                                            FUN = fun.data,
                                            na.rm = na.rm,
                                            na.action = na.pass)
                                summary
                                }
                             )

update_r <-
  function(r,
           data,
           interp_limit) {
    abs_x_max <- max(abs(data$x),
                     na.rm = TRUE)
    abs_y_max <- max(abs(data$y),
                     na.rm = TRUE)
    r <- switch(interp_limit,
                 "head" = sqrt(abs_x_max^2 + abs_y_max^2),
                 "skirt" = r) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
    r
  }
