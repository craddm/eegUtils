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
                                            method) {

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
                      r = 95,
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
                                                      r) {
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

#' @noRd
biharmonic <- function(data,
                       grid_res,
                       interp_limit,
                       r = 95) {

  # max_dim <- max(c(abs(data$x),
  #                  abs(data$y)))
  abs_y_max <- max(abs(data$y),
                   na.rm = TRUE)
  abs_x_max <- max(abs(data$x),
                   na.rm = TRUE)

  x_min <- min(data$x, na.rm = TRUE) * 2.5
  x_max <- max(data$x, na.rm = TRUE) * 2.5
  y_min <- min(data$y, na.rm = TRUE) * 2.5
  y_max <- max(data$y, na.rm = TRUE) * 2.5

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
  g <- (d^2) * (log(d) - 1) #Green's function
  diag(g) <- 0
  weights <- qr.solve(g, data$fill)
  xy <- t(xy)
  # Remind me to make this code readable at some point.

  outmat <- numeric(grid_res^2)

  one_fun <- function(xo, yo) {

    tmp <-
      matrix(complex(real = xo,
                     imaginary = yo),
             nrow = grid_res,
             ncol = grid_res)
    tmp_x <- numeric(length(xy))

    for (i in seq(length(tmp))) {
      tmp_x <- (abs(tmp[i] - xy)^2) * (log(abs(tmp[i] - xy)) - 1)
      tmp_x <- ifelse(is.nan(tmp_x),
                      0,
                      tmp_x)
      outmat[i] <- tmp_x %*% weights
    }
    outmat
  }

  outmat <- one_fun(xo, yo)

  dim(outmat) <- c(grid_res,
                   grid_res)
  data <- data.frame(x = xo[, 1],
                     outmat)
  names(data)[1:length(yo[1, ]) + 1] <- yo[1, ]

  data <-
    tidyr::pivot_longer(data,
                        cols = !x,
                        names_to = "y",
                        values_to = "fill",
                        names_transform = list(y = as.numeric))

  if (identical(interp_limit,
                "head")) {
    circ_scale <- sqrt(abs_x_max^2 + abs_y_max^2) * 1.01
    } else {
     circ_scale <- sqrt(
       ((abs_x_max/2)^2 + abs_x_max^2) +
         ((abs_y_max/2)^2 + abs_y_max^2)
       )
  }
  data[sqrt(data$x ^ 2 + data$y ^ 2) <= circ_scale, ]

}


#' @noRd
fit_gam_topo <- function(data,
                         grid_res,
                         interp_limit) {

  abs_y_max <- max(abs(data$y),
                   na.rm = TRUE)
  y_max <- max(data$y,
               na.rm = TRUE) * 2.5

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


  if (identical(interp_limit,
                "head")) {
    circ_scale <- abs_y_max * 1.1
  } else {
    circ_scale <- y_max / 1.8
  }

  data$incircle <- sqrt(data$x ^ 2 + data$y ^ 2) < circ_scale
  data[data$incircle, ]
}

#' Get an interpolated grid
#'
#' @param data Object to be interpolated
#' @param ... Other arguments
#' @export
get_scalpmap <- function(data,
                         ...) {
  UseMethod("get_scalpmap", data)
}

#' @export
get_scalpmap.default <- function(data,
                                 ...) {
  stop("Not implemented for objects of class ", class(data))
}

#' @param data An object of class 'data.frame`
#' @param method "biharmonic" or "gam"
#' @param grid_res grid resolution
#' @param interp_limit interpolate up to the "head" or add a "skirt"
#' @param quantity "amplitude"
#' @param facets Any facets you plan to use
#' @describeIn get_scalpmap data.frame
#' @export
get_scalpmap.data.frame <- function(data,
                                    method = "biharmonic",
                                    grid_res = 100,
                                    interp_limit = "skirt",
                                    quantity = "amplitude",
                                    facets = NULL,
                                    ...) {

  facets <- rlang::enexpr(facets)
  method <- match.arg(method,
                      c("biharmonic",
                        "gam",
                        "Biharmonic"))
  method <- switch(method,
                   "Biharmonic" = "biharmonic",
                   method)

  if (!is.null(facets)) {
    tmp <- dplyr::group_by(tmp,
                           {{facets}})
  } else {
    tmp <- data
  }

  if (!is.null(facets)) {
    tmp <- dplyr::group_nest(tmp, {{facets}})
    tmp <- dplyr::mutate(tmp,
                         topos = map(data,
                                     ~biharmonic(.,
                                                 grid_res = grid_res,
                                                 interp_limit = interp_limit)),
    )
    tmp <- dplyr::select(tmp, -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit)
      )
  }
  smooth

}

#' @export
get_scalpmap.eeg_epochs <- function(data,
                                    method = "biharmonic",
                                    grid_res = 100,
                                    interp_limit = "skirt",
                                    quantity = "amplitude",
                                    facets = NULL,
                                     ...) {

  facets <- rlang::enexpr(facets)
  # only useful if passed eeg objects
  check_locs <- no_loc_chans(channels(data))

  if (!is.null(check_locs)) {
    data <- select(data, -check_locs)
  }

  tmp <- as.data.frame(data)
  if (!is.null(facets)) {
    tmp <- dplyr::group_by(tmp, {{facets}})
  }
  #tmp <- group_by(tmp, event_label)
  tmp <- dplyr::summarise_at(tmp,
                             channel_names(data),
                             .funs = mean)
  tmp <- tidyr::gather(tmp,
                       electrode,
                       amplitude,
                       -{{facets}})
  tmp <- dplyr::left_join(tmp,
                          channels(data),
                          by = "electrode")
  tmp <- dplyr::rename(tmp, fill = {{quantity}})

  if (!is.null(facets)) {
    tmp <- dplyr::group_nest(tmp, {{facets}})
    tmp <- dplyr::mutate(tmp,
                         topos = map(data,
                                     ~biharmonic(.,
                                                 grid_res = grid_res,
                                                 interp_limit = interp_limit)),
                         )
    tmp <- dplyr::select(tmp, -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit)
            )
  }
  smooth
}

#' @export
get_scalpmap.eeg_ICA <- function(data,
                                 method = "biharmonic",
                                 grid_res = 100,
                                 interp_limit = "skirt",
                                 quantity = "amplitude",
                                 facets = component,
                                 verbose = FALSE,
                                 ...) {

  facets <- rlang::enexpr(facets)

  tmp <- as.data.frame(data,
                       mixing = TRUE,
                       long = TRUE)

  if (any(is.na(tmp$x))) {
    tmp <- tmp[!is.na(tmp$x), ]
    if (verbose) {
      warning("Removing channels with no location.")
    }
  }

  tmp <- dplyr::rename(tmp,
                       fill = {{quantity}})

  if (!is.null(facets)) {

    tmp <- dplyr::group_nest(tmp,
                             component)
    tmp <- dplyr::mutate(tmp,
                         topos = map(data,
                                     ~biharmonic(.,
                                                 grid_res = grid_res,
                                                 interp_limit = interp_limit)),
    )
    tmp <- dplyr::select(tmp,
                         -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit)
      )
  }
  smooth

}

#' @noRd
no_loc_chans <- function(chaninfo) {

  if (any(is.na(chaninfo$x))) {
    return(chaninfo$electrode[is.na(chaninfo$x)])
  }
  NULL
}


stat_summary_by_fill <- function(mapping = NULL,
                                 data = NULL,
                                 geom = "raster",
                                 position = "identity",
                                  ...,
                                 fun.data = NULL,
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = TRUE) {
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
