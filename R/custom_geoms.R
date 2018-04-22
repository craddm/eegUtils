
#' StatBiharmonic
#'
#'
#'

StatBiharmonic <- ggproto("StatBiharmonic", Stat,
                          required_aes = c("x", "y", "fill"),

                          compute_group = function(data, scales) {
                            data <- aggregate(fill ~ x + y, data = data, FUN = mean)

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
stat_biharmonic <- function(mapping = NULL, data = NULL, geom = "raster",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatBiharmonic, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



StatHead <- ggplot2::ggproto("StatHead",
                             Stat,
                             default_aes = aes(colour = ..id..),
                             compute_group = function(data, scales) {
                               r <- 1
                               head_shape <- data.frame(x = r * cos(seq(0, 2 * pi, length.out = 101)),
                                                        y = r * sin(seq(0, 2 * pi, length.out = 101)),
                                                        id = "head")

                               nose <- data.frame(x = c(head_shape$x[[23]],
                                                        head_shape$x[[26]],
                                                        head_shape$x[[29]]),
                                                  y = c(head_shape$y[[23]],
                                                        head_shape$y[[26]] * 1.1,
                                                        head_shape$y[[29]]),
                                                  id = "nose")

                               ears <- data.frame(x = c(head_shape$x[[4]],
                                                        head_shape$x[[97]],
                                                        head_shape$x[[48]],
                                                        head_shape$x[[55]]),
                                                  y = c(head_shape$y[[4]],
                                                        head_shape$y[[97]],
                                                        head_shape$y[[48]],
                                                        head_shape$y[[55]]),
                                                  id = "ears")
                               data <- rbind(head_shape, nose, ears)
                               data
                             }
)

#' @inheritParams geom_path
stat_head <- function(mapping = NULL, data = NULL, geom = "path",
                      position = "identity", na.rm = FALSE,
                      show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatHead, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}




#' StatScalpmap
#'

StatScalpmap <- ggproto("StatScalpmap", Stat,
                          required_aes = c("x", "y", "fill"),

                          compute_group = function(data, scales) {
                            data <- aggregate(fill ~ x + y, data = data, FUN = mean)

                            x_min <- floor(min(data$x))
                            x_max <- ceiling(max(data$x))
                            y_min <- floor(min(data$y))
                            y_max <- ceiling(max(data$y))

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
                            data <- data.frame(x = xo[, 1], outmat)
                            names(data)[1:length(yo[1, ]) + 1] <- yo[1, ]
                            data <- tidyr::gather(data,
                                                  key = y,
                                                  value = fill,
                                                  -x,
                                                  convert = TRUE)
                            radius <- max(sqrt(data$x ^2 + data$y ^2), na.rm = T) * .57
                            data$incircle <- sqrt(data$x ^ 2 + data$y ^ 2) < radius

                            data[!data$incircle, ] <- NA
                            data
                          }
)

#' @inheritParams geom_raster
stat_scalpmap <- function(mapping = NULL, data = NULL, geom = "raster",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatScalpmap, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
