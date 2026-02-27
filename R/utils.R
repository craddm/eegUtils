#' Convert to 3d matrix
#'
#' @param data data to be converted
#' @param ... additional parameters
#' @keywords internal
conv_to_mat <- function(data, ...) {
  UseMethod("conv_to_mat", data)
}

#' @describeIn conv_to_mat Default method
conv_to_mat.default <- function(data, ...) {
  stop("Not implemented for objects of class", class(data))
}

#' @describeIn conv_to_mat Convert `eeg_epochs` to 3D matrix
conv_to_mat.eeg_epochs <- function(data, ...) {
  n_epochs <- length(unique(data$timings$epoch))
  n_channels <- ncol(data$signals)
  n_times <- length(unique(data$timings$time))
  data <- array(as.matrix(data$signals),
                dim = c(n_times, n_epochs, n_channels))
  data
}

is.true <- function(x) {
  !is.na(x) & x
}

#' Reflect-pad a vector
#'
#' Pads a vector by reflecting the signal at both boundaries. The reflection is
#' odd (anti-symmetric): `2*x[1] - rev(x[2:(n+1)])` at the start and
#' `2*x[end] - rev(x[(end-n):(end-1)])` at the end. This matches the padding
#' used by SciPy's `filtfilt` and avoids the edge discontinuities caused by
#' zero-padding.
#'
#' @param x Numeric vector to pad.
#' @param n Number of samples to add at each end.
#' @return A numeric vector of length `length(x) + 2*n`.
#' @keywords internal
reflect_pad <- function(x, n) {
  len <- length(x)
  if (n >= len) {
    stop("Pad length must be less than signal length.")
  }
  start_pad <- 2 * x[1] - rev(x[2:(n + 1)])
  end_pad <- 2 * x[len] - rev(x[(len - n):(len - 1)])
  c(start_pad, x, end_pad)
}

#' Pad a vector with zeros
#'
#' Pads a vector with zeros at the beginning and end.
#'
#' @param x vector to pad
#' @param n number of zeros to pad
#' @param startval value to add at the start of the vector - defaults to zero
#' @param endval value to add at the end of the vector - defaults to zero
#' @keywords internal
pad <- function(x,
                n,
                startval = 0,
                endval = 0) {
  c(rep(startval, n),
    x,
    rep(endval, n))
  }

#' Unpad a vector
#'
#' remove a specified number of points from the beginning and end of a vector
#' @param x Vector to unpad/trim
#' @param n Number of points to remove
#' @keywords internal
unpad <- function(x, n) {
  start <- n + 1
  end <- length(x) - n
  x <- x[start:end]
  x
}

#' Get number of samples
#' @param .data Object to get total number of sampling points from
#' @keywords internal
samples <- function(.data) {
  nrow(.data$signals)
}

#' Generate circle as radians
#' @keywords internal
circ_rad_fun <- function() {
  seq(0,
      2 * pi,
      length.out = 101)
}

#' Update radius
#' @keywords internal
update_r <-
  function(r = 95,
           data,
           interp_limit) {

    max_elec <- calc_max_elec(data)
    r <- switch(interp_limit,
                "head" = min(max_elec * 1.10, max_elec + 10),
                "skirt" = r,
                "convex_hull" = min(max_elec * 1.10, max_elec + 10)) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
    r
  }

#' Calculate maximum electrode distance from origin.
#' @keywords internal
calc_max_elec <- function(data) max(sqrt(data$x^2 + data$y^2), na.rm = TRUE)

#' Test whether points lie inside the convex hull of electrode positions
#'
#' Uses the cross-product winding method for convex polygons. Optionally
#' expands the hull outward by a small buffer so that electrode positions
#' themselves are not clipped.
#'
#' @param elec_x,elec_y Numeric vectors of electrode coordinates.
#' @param grid_x,grid_y Numeric vectors of grid coordinates to test.
#' @param buffer Fractional expansion of the hull (default 0.10, i.e. 10 percent).
#' @return A logical vector the same length as `grid_x`.
#' @keywords internal
point_in_hull <- function(elec_x, elec_y, grid_x, grid_y, buffer = 0.10) {

  hull_idx <- grDevices::chull(elec_x, elec_y)
  hx <- elec_x[hull_idx]
  hy <- elec_y[hull_idx]

  # expand hull outward from centroid
  cx <- mean(hx)
  cy <- mean(hy)
  hx <- cx + (1 + buffer) * (hx - cx)
  hy <- cy + (1 + buffer) * (hy - cy)

  n <- length(hx)
  inside <- rep(TRUE, length(grid_x))

  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    # cross product of edge vector with point-to-vertex vector
    cross <- (hx[j] - hx[i]) * (grid_y - hy[i]) -
             (hy[j] - hy[i]) * (grid_x - hx[i])
    inside <- inside & (cross <= 0)
  }
  inside
}

