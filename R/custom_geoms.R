# GeomCircle <- ggplot2::ggproto("GeomCircle", ggplot2::Geom,
#   required_aes = c("x", "y", "radius"),
#   default_aes = ggplot2::aes(colour = "black", fill=NA,
#                              alpha=NA, linewidth=1, linetype="solid"))
#
#
# geom_circle <- function(mapping = NULL, data = NULL, stat = "circle",
#                         position = "identity", n = 360, expand = 0, radius = 0, na.rm = FALSE,
#                         show.legend = NA, inherit.aes = TRUE, ...) {
#   ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomCircle,
#         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#         params = list(n = n, na.rm = na.rm, ...))
#
# }


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

                            #xo <- seq(x_min + x_min / 4, x_max + x_max / 4, length = 80)
                            #yo <- seq(y_min + y_min / 4, y_max + y_max / 4, length = 80)
                            xo <- seq(x_min, x_max, length = 80)
                            yo <- seq(y_min, y_max, length = 80)

                            xo <- matrix(rep(xo, 80),
                                         nrow = 80,
                                         ncol = 80)

                            yo <- t(matrix(rep(yo, 80),
                                           nrow = 80,
                                           ncol = 80))

                            #max_dim <- max(abs(data$x), abs(data$y))
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

