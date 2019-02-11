#' @method fortify eeg_epochs
fortify.eeg_epochs <- function(model,
                               data,
                               ...) {
  tibble::as_tibble(as.data.frame(model,
                                  long = TRUE,
                                  stringsAsFactors = FALSE))
}


fortify.eeg_data <- function(model,
                               data,
                               ...) {
  as.data.frame(model,
                long = TRUE ,
                stringsAsFactors = FALSE)
}

fortify.eeg_ICA <- function(model,
                            data,
                            ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}


fortify.eeg_tfr <- function(model,
                            data,
                            ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}

fortify.eeg_evoked <- function(model,
                               data,
                               ...) {
  as.data.frame(model,
                long = TRUE,
                stringsAsFactors = FALSE)
}


#' StatBiharmonic
#' @noRd

StatBiharmonic <- ggproto("StatBiharmonic",
                          Stat,
                          required_aes = c("x", "y", "fill"),

                          compute_group = function(data,
                                                   scales) {
                            data <- aggregate(fill ~ x + y,
                                              data = data,
                                              FUN = mean)

                            x_min <- min(data$x)
                            x_max <- max(data$x)
                            y_min <- min(data$y)
                            y_max <- max(data$y)

                            xo <- seq(x_min, x_max, length = 80)
                            yo <- seq(y_min, y_max, length = 80)

                            xo <- matrix(rep(xo, 80),
                                         nrow = 80,
                                         ncol = 80)

                            yo <- t(matrix(rep(yo, 80),
                                           nrow = 80,
                                           ncol = 80))

                            xy_coords <- unique(data[, c("x", "y")])

                            xy <- xy_coords[, 1] + xy_coords[, 2] * sqrt(as.complex(-1))

                            d <- matrix(rep(xy, length(xy)),
                                        nrow = length(xy),
                                        ncol = length(xy))

                            d <- abs(d - t(d))
                            diag(d) <- 1
                            g <- (d ^ 2) * (log(d) - 1) #Green's function
                            diag(g) <- 0
                            weights <- qr.solve(g, data$fill)
                            xy <- t(xy)

                            outmat <-
                              purrr::map(xo + sqrt(as.complex(-1)) * yo,
                                         function (x) (abs(x - xy) ^ 2) *
                                           (log(abs(x - xy)) - 1) ) %>%
                              rapply(function (x) ifelse(is.nan(x), 0, x), how = "replace") %>%
                              purrr::map_dbl(function (x) x %*% weights)

                            dim(outmat) <- c(80, 80)
                            out_df <- data.frame(x = xo[, 1], outmat)
                            names(out_df)[1:length(yo[1, ]) + 1] <- yo[1, ]
                            data <- tidyr::gather(out_df,
                                                  key = y,
                                                  value = fill,
                                                  -x,
                                                  convert = TRUE)
                            data
                          }
)

#' @inheritParams geom_raster
stat_biharmonic <- function(mapping = NULL,
                            data = NULL,
                            geom = "raster",
                            position = "identity",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE,
                            ...) {
  ggplot2::layer(
    stat = StatBiharmonic,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



#' StatScalpmap
#' @noRd

StatScalpmap <- ggplot2::ggproto("StatScalpmap",
                                 Stat,
                                 required_aes = c("x",
                                                  "y",
                                                  "fill"),

                        compute_group = function(data,
                                                 scales) {

                           data <- aggregate(fill ~ x + y,
                                             data = data,
                                             FUN = mean)

                           x_min <- floor(min(data$x)) * 2
                           x_max <- ceiling(max(data$x)) * 2
                           y_min <- floor(min(data$y)) * 2
                           y_max <- ceiling(max(data$y)) * 2

                           xo <- seq(x_min,
                                     x_max,
                                     length = 80)
                           yo <- seq(y_min,
                                     y_max,
                                     length = 80)
                           xo <- matrix(rep(xo,
                                            80),
                                        nrow = 80,
                                        ncol = 80)

                           yo <- t(matrix(rep(yo, 80),
                                          nrow = 80,
                                          ncol = 80))

                          xy_coords <- unique(data[, c("x", "y")])

                          xy <- xy_coords[, 1] + xy_coords[, 2] * sqrt(as.complex(-1))

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

                          dim(outmat) <- c(80, 80)
                          data <- data.frame(x = xo[, 1],
                                             outmat)
                          names(data)[1:length(yo[1, ]) + 1] <- yo[1, ]
                          data <- tidyr::gather(data,
                                                key = y,
                                                value = fill,
                                                -x,
                                                convert = TRUE)

                          circ_scale <- y_max / 1.5
                          data$incircle <- sqrt(data$x ^ 2 + data$y ^ 2) < circ_scale
                          data[data$incircle, ]
                        }
)

#' @inheritParams geom_raster
stat_scalpmap <- function(mapping = NULL,
                          data = NULL,
                          geom = "raster",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE,
                          ...) {
  ggplot2::layer(
    stat = StatScalpmap,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  ...)
  )
}

geom_topo <- function(mapping = NULL,
                      data = NULL,
                      stat = "identity",
                      position = "identity",
                      show.legend = NA,
                      na.rm = TRUE,
                      inherit.aes = TRUE,
                      interpolate = TRUE,
                      chan_markers = "point",
                      fill = NA,
                      size = NA,
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
                                    ...)
                      ),
       ggplot2::layer(geom = GeomPath,
                      stat = StatMask,
                      data = data,
                      mapping = mapping,
                      position = position,
                      inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm,
                                    colour = "white",
                                    size = rel(6.5),
                                    ...)),
       if (chan_markers == "point") {
         ggplot2::layer(geom = GeomPoint,
                        data = data,
                        stat = stat,
                        mapping = mapping,
                        position = position,
                        show.legend = show.legend,
                        inherit.aes = inherit.aes,
                        params = list(na.rm = na.rm,
                                      fill = fill,
                                      size = size,
                                      ...))}
  )
}

GeomTopo <- ggplot2::ggproto("GeomTopo",
                             GeomRaster)

StatHead <- ggplot2::ggproto("StatHead",
                             Stat,
                             compute_group = function(data,
                                                      scales) {

                               y_lim <- max(data$y) * 1.1
                               make_head(r = y_lim)
                             }
                             )

GeomHead <- ggplot2::ggproto("GeomHead",
                              GeomPath)

geom_head <- function(mapping = NULL,
                      data = NULL,
                      show.legend = NA,
                      na.rm = TRUE,
                      inherit.aes = TRUE,
                      ...) {

  ggplot2::layer(geom = GeomHead,
                 data = data,
                 mapping = mapping,
                 stat = StatHead,
                 position = PositionIdentity,
                 inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm,
                               ...))
}



# stat_mask <- function(mapping = NULL,
#                       data = NULL,
#                       stat = "identity",
#                       position = "identity",
#                       show.legend = NA,
#                       na.rm = TRUE,
#                       inherit.aes = TRUE,
#                       geom = "mask",
#                       ...) {
#
#   ggplot2::layer(stat = StatMask,
#                  data = data,
#                  mapping = mapping,
#                  geom = geom,
#                  position = position,
#                  show.legend = show.legend,
#                  inherit.aes = inherit.aes,
#                  params = list(na.rm = na.rm,
#                                ...)
#                  )
#
# }

StatMask <-
  ggplot2::ggproto("StatMask",
                   Stat,
                   compute_group = function(data,
                                            scales) {

                     scale_fac <- 1.4 * max(data$y)
                     data <- data.frame(x = scale_fac * cos(circ_rad_fun()),
                                        y = scale_fac * sin(circ_rad_fun()))
                     data

                   }
)

geom_mask <- function(mapping = NULL,
                      data = NULL,
                      stat = "identity",
                      position = "identity",
                      show.legend = NA,
                      na.rm = FALSE,
                      colour = "white",
                      size = rel(6.5),
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
                               ...))
}

GeomEars <- ggplot2::ggproto("GeomEars",
                             GeomCurve)

StatEars <- ggplot2::ggproto("StatEars",
                             Stat,
                             compute_group = function(data, scales) {

                             })

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

  # head_out <- list("head_shape" = head_shape,
  #                  "nose" = nose,
  #                  "ears" = ears)
  head_out <- rbind(head_shape, nose)
  head_out
}

make_l_ear <- function(r) {

  head_shape <- data.frame(x = r * cos(circ_rad_fun()),
                           y = r * sin(circ_rad_fun()))
  left_ear <- data.frame(x = head_shape$x[[4]],
                         xend = head_shape$x[[97]],
                         y = head_shape$y[[4]],
                         yend = head_shape$y[[97]])
  left_ear
}

make_r_ear <- function(r) {
  right_ear <- data.frame(x = c(head_shape$x[[48]],
                                head_shape$x[[55]]),
                          y = c(head_shape$y[[48]],
                                head_shape$y[[55]]))
}
