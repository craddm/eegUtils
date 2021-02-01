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
#' @param method "biharmonic" or "gam" smooth
#' @param grid_res Grid resolution
#' @param interp_limit Interpolate up to the "head" or add a "skirt"
#' @param quantity Quantity to be plotted. Defaults to "amplitude".
#' @param facets Any facets you plan to use
#' @param r Size of headshape
#' @describeIn get_scalpmap data.frame
#' @export
get_scalpmap.data.frame <- function(data,
                                    method = "biharmonic",
                                    grid_res = 100,
                                    interp_limit = "skirt",
                                    quantity = "amplitude",
                                    facets = NULL,
                                    r = 95,
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
                                                 interp_limit = interp_limit,
                                                 r = r)),
    )
    tmp <- dplyr::select(tmp, -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r)
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
                                    r = 95,
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
                                                 interp_limit = interp_limit,
                                                 r = r)),
    )
    tmp <- dplyr::select(tmp, -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r)
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
                                 r = 95,
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
                                                 interp_limit = interp_limit,
                                                 r = r)),
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
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r)
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


#' @noRd
biharmonic <- function(data,
                       grid_res,
                       interp_limit,
                       r = NULL) {

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

    if (is.null(r)) {
      circ_scale <- sqrt(abs_x_max^2 + abs_y_max^2) * 1.01
    } else {
      circ_scale <- r * 1.01
    }

  } else {
    max_elec <- max(abs_x_max,
                    abs_y_max,
                    na.rm = TRUE)
    if (max_elec > r) {
      circ_scale <- (max_elec - r) + max_elec
    } else {
      circ_scale <- r * 1.25
    }
  }
  data[sqrt(data$x ^ 2 + data$y ^ 2) <= circ_scale, ]

}


#' @noRd
fit_gam_topo <- function(data,
                         grid_res,
                         interp_limit,
                         r) {

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
