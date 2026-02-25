#' Calculate an interpolated scalpmap
#'
#' Calculate and return an interpolated scalpmap from a `data.frame` or
#' `eeg_epochs` object. This provides the basis for a raster plot, as used to
#' create topographical plots.
#'
#' @param data Object to be interpolated
#' @param ... Other arguments
#' @return A `tibble` with the columns `x`, `y`, `fill`
#' @examples
#' head(get_scalpmap(demo_epochs))
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
#' @param interp_limit Interpolate up to the "head", add a "skirt", or clip to
#'   the "convex_hull" of electrode positions
#' @param quantity Quantity to be plotted. Defaults to "amplitude".
#' @param facets Any facets you plan to use
#' @param r Size of headshape
#' @param k Degrees of freedom used for spline when using `method = gam`.
#'   Defaults to -1, which attempts to automatically determine k.
#' @describeIn get_scalpmap Get scalpmap from a `data.frame`
#' @export
get_scalpmap.data.frame <- function(data,
                                    method = "biharmonic",
                                    grid_res = 100,
                                    interp_limit = "skirt",
                                    quantity = "amplitude",
                                    facets = NULL,
                                    r = 95,
                                    k = -1,
                                    ...) {

  method <- match.arg(tolower(method), c("biharmonic", "gam"))

  if (!is.null(facets)) {
  #if (!is.null(rlang::enquo(facets))) {
    facets <- rlang::enexpr(facets)
    tmp <- dplyr::group_by(data,
                           dplyr::across({{facets}}))
    tmp <- dplyr::group_nest(tmp)
    tmp <-
      dplyr::mutate(tmp,
                    topos = map(data,
                                ~switch(method,
                                        biharmonic(.,
                                                   grid_res = grid_res,
                                                   interp_limit = interp_limit,
                                                   r = r),
                                        gam = fit_gam_topo(.,
                                                           grid_res = grid_res,
                                                           interp_limit = interp_limit,
                                                           r = r,
                                                           k = k)
                                        )
                                )
                    )
    tmp <- dplyr::select(tmp, -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    tmp <- data
    smooth <-
      switch(method,
             biharmonic = biharmonic(tmp,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r),
             gam = fit_gam_topo(tmp,
                                grid_res = grid_res,
                                interp_limit = interp_limit,
                                r = r,
                                k = k)
      )
  }
  smooth

}

#' @describeIn get_scalpmap Get scalpmap from an `eeg_epochs` object.
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

  calc_map(data = data,
           method = method,
           grid_res = grid_res,
           interp_limit = interp_limit,
           quantity = quantity,
           facets = facets,
           r = r)
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
  tmp <- tmp[stats::complete.cases(tmp),]
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

calc_map <- function(data,
                     method,
                     grid_res,
                     interp_limit,
                     quantity,
                     facets,
                     r) {

  tmp <- as.data.frame(data)

  if (!is.null(facets)) {
    tmp <- dplyr::group_by(tmp,
                           dplyr::across(!!{{ facets }}))
  }
  tmp <- dplyr::summarise(tmp,
                          dplyr::across(channel_names(data),
                                        .fns = mean),
                          .groups = "keep")
  tmp <- tidyr::pivot_longer(tmp,
                             cols = channel_names(data),
                             names_to = "electrode",
                             values_to = "fill")
  tmp <- dplyr::left_join(tmp,
                          subset(channels(data),
                                 select = c(electrode, x, y)))
  if (!is.null(facets)) {
    tmp <- dplyr::group_nest(tmp)
    if (identical(method, "biharmonic")) {
      tmp <-
        dplyr::mutate(tmp,
                      topos = map(data,
                                  ~biharmonic(.,
                                              grid_res = grid_res,
                                              interp_limit = interp_limit,
                                              r = r)))
    } else {
      tmp$topos <- map(tmp$data,
                       ~fit_gam_topo(.,
                                     grid_res = grid_res,
                                     interp_limit = interp_limit,
                                     r = r)
                       )
    }

    tmp <- subset(tmp,
                  select = -data)
    smooth <- tidyr::unnest(tmp,
                            cols = topos)
  } else {
    tmp <- subset(tmp,
                  select = c(x, y, electrode, fill))
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


#' Clip interpolated grid data to the desired boundary
#' @noRd
clip_to_limit <- function(data, interp_limit, r, max_elec, elec_x, elec_y) {

  if (identical(interp_limit, "convex_hull")) {
    in_hull <- point_in_hull(elec_x, elec_y, data$x, data$y)
    return(data[in_hull, ])
  }

  circ_scale <- if (identical(interp_limit, "head")) {
    if (is.null(r)) max_elec * 1.01 else r * 1.01
  } else {
    # skirt: 20% or 20 mm buffer past furthest electrode, whichever is smaller
    ref <- max(r, max_elec)
    min(ref * 1.20, ref + 20)
  }

  data[sqrt(data$x^2 + data$y^2) <= circ_scale, ]
}

#' @noRd
biharmonic <- function(data,
                       grid_res,
                       interp_limit,
                       r = NULL) {

  max_elec <- calc_max_elec(data)

  x_min <- r * -2
  x_max <- r * 2
  y_min <- r * -2
  y_max <- r * 2

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

  # Vectorised Green's function evaluation over the full grid
  tmp <- matrix(complex(real = xo,
                        imaginary = yo),
                nrow = grid_res,
                ncol = grid_res)
  d_grid <- abs(outer(as.vector(tmp), as.vector(xy), `-`))
  g_grid <- d_grid^2 * (log(d_grid) - 1)
  g_grid[is.nan(g_grid)] <- 0
  outmat <- matrix(drop(g_grid %*% weights),
                   nrow = grid_res,
                   ncol = grid_res)

  data <- data.frame(x = xo[, 1],
                     outmat)
  names(data)[1:length(yo[1, ]) + 1] <- yo[1, ]

  data <-
    tidyr::pivot_longer(data,
                        cols = !x,
                        names_to = "y",
                        values_to = "fill",
                        names_transform = list(y = as.numeric))

  clip_to_limit(data, interp_limit, r, max_elec, xy_coords$x, xy_coords$y)

}


#' @noRd
fit_gam_topo <- function(data,
                         grid_res,
                         interp_limit,
                         r,
                         k = -1) {

  max_elec <- calc_max_elec(data)
  elec_xy <- unique(data[, c("x", "y")])

  spline_smooth <- mgcv::gam(fill ~ s(x,
                                      y,
                                      bs = "tp",
                                      k = k),
                             data = data)

  data <- expand.grid(x = seq(max_elec * -1.5,
                              max_elec * 1.5,
                              length = grid_res),
                      y = seq(max_elec * -1.5,
                              max_elec * 1.5,
                              length = grid_res))
  data$fill <- stats::predict(spline_smooth,
                              data,
                              type = "response")

  clip_to_limit(data, interp_limit, r, max_elec, elec_xy$x, elec_xy$y)
}
