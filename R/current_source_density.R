#' Channel interpolation
#'
#' Interpolate EEG channels using a spherical spline (Perrin et al., 1989). The
#' data must have channel locations attached.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param bad_elecs Name(s) of electrode(s) to interpolate.
#' @param ... Other parameters passed to the functions.
#' @export

interp_elecs <- function(data, bad_elecs, ...) {
  UseMethod("interp_elecs", data)
}

#' @describeIn interp_elecs Interpolate EEG channel(s)
#' @export
interp_elecs.eeg_data <- function(data,
                                  bad_elecs, ...) {

  if (is.null(data$chan_info)) {
    stop("No channel locations found.")
  }

  if (all(c("theta", "radius") %in% names(data$chan_info))) {
    xyz_coords <- pol_to_sph(data$chan_info$theta,
                             data$chan_info$radius)
  } else {
    xyz_coords <- pol_to_sph(data$chan_info$pol_theta,
                             data$chan_info$pol_radius)
  }
  rads <- sqrt(xyz_coords$x ^ 2 + xyz_coords$y ^ 2 + xyz_coords$z ^ 2)
  xyz_coords <- xyz_coords / rads

  bad_select <- toupper(data$chan_info$electrode) %in% toupper(bad_elecs)
  xyz_bad <- xyz_coords[bad_select, ]
  xyz_good <- xyz_coords[!bad_select, ]

  if (nrow(xyz_bad) == 0) {
    return(data)
  }

  sigs_select <- toupper(names(data$signals)) %in%
    toupper(data$chan_info$electrode)

  sigs <- data$signals[sigs_select]

  bad_cols <- toupper(names(sigs)) %in% toupper(bad_elecs)

  mm <- spheric_spline(xyz_good, xyz_bad, sigs[!bad_cols])

  if (ncol(mm) == 1) {
    mm <- as.vector(mm)
  }
  sigs[bad_cols] <- mm

  data$signals[, sigs_select] <- sigs
  data
}


#' Calculate a spherical spline smooth for interpolation of electrodes
#'
#' @author Matt Craddock \email{matt@mattcraddock.com}
#'
#' @param good_elecs Electrodes with positions that do not need interpolation
#' @param bad_elecs Electrodes to be interpolated
#' @param data Raw data
#' @noRd

spheric_spline <- function(good_elecs, bad_elecs, data) {

  lec <- compute_g(good_elecs, good_elecs)
  ph <- compute_g(bad_elecs, good_elecs)

  meandata <- rowMeans(data)
  data <- data - meandata
  data$CC <- 0

  C <- as.matrix(data) %*% t(MASS::ginv(rbind(lec, 1), tol = 0))

  mmm <- matrix(NA, nrow = nrow(bad_elecs), ncol = nrow(data))

  if (nrow(bad_elecs) > 1) {
    for (i in seq(1, nrow(bad_elecs))) {
      mmm[i, ] <- rowSums(C * t(matrix(rep(ph[i, ],
                                           nrow(data)),
                                       nrow = ncol(C))))
    }
  } else {
    mmm <- rowSums(C * t(matrix(rep(ph, nrow(data)), nrow = ncol(C))))
  }

  mmm <- mmm + matrix(rep(meandata, nrow(bad_elecs)), nrow = nrow(bad_elecs))
  t(mmm) # output should be transposed to be in appropriate columns

}

#' Compute the g function for two sets of locations of channel locations on the
#' unit sphere.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param xyz_coords A set of electrode locations on a unit sphere.
#' @param xyz_elecs A set of electrode locations on a unit sphere.
#' @param m Interpolation constant (higher = less flexible)
#' @importFrom pracma legendre
#' @noRd

compute_g <- function(xyz_coords,
                      xyz_elecs,
                      m = 4) {

  EI <- 1 - sqrt((matrix(rep(xyz_coords[, 1],
                             nrow(xyz_elecs)),
                         nrow = nrow(xyz_elecs)) -
                    matrix(rep(xyz_elecs[, 1],
                               each = nrow(xyz_coords)),
                           ncol = nrow(xyz_coords))) ^ 2 +
                   (matrix(rep(xyz_coords[, 2],
                               nrow(xyz_elecs)),
                           nrow = nrow(xyz_elecs)) -
                      matrix(rep(xyz_elecs[, 2],
                                 each = nrow(xyz_coords)),
                             ncol = nrow(xyz_coords))) ^ 2 +
                   (matrix(rep(xyz_coords[, 3],
                               nrow(xyz_elecs)),
                           nrow = nrow(xyz_elecs)) -
                      matrix(rep(xyz_elecs[, 3],
                                 each = nrow(xyz_coords)),
                             ncol = nrow(xyz_coords))) ^ 2)

  dim(EI) <- c(nrow(xyz_coords), nrow(xyz_elecs))

  g <- 0

  for (i in seq(1, 7)) {
    poly_xy <- pracma::legendre(i, EI)
    dim(poly_xy) <- c(i + 1,
                      nrow(EI),
                      ncol(EI))
    g <- g + ((2 * i + 1) / (i ^ m * (i + 1) ^ m)) * poly_xy[1, , ]
  }
  g <- g / (4 * pi)
}




#' Calculate current source density
#'
#' @param data eeg_data object
#' @param m smoothing constraint (higher = more rigid)
#' @param smoothing lambda constant
#' @param scaling Unsure - don't think my coords are typically on unit sphere...
#' @noRd

compute_csd <- function(data, m = 4, smoothing = 1e-05, scaling = 10) {

  orig_elecs <- names(data$signals)
  data_chans <- toupper(names(data$signals)) %in% toupper(data$chan_info$electrode)

  if (any(!data_chans)) {
    data <- reref_eeg(data,
                      exclude = names(data$signals)[!data_chans])
  } else {
    data <- reref_eeg(data)
  }
  if (all(c("theta", "radius") %in% names(data$chan_info))) {
    xyz_coords <- pol_to_sph(data$chan_info$theta,
                             data$chan_info$radius)
  } else {
    xyz_coords <- pol_to_sph(data$chan_info$pol_theta,
                             data$chan_info$pol_radius)
  }
  rads <- sqrt(xyz_coords$x ^ 2 + xyz_coords$y ^ 2 + xyz_coords$z ^ 2)
  xyz_coords <- xyz_coords / rads
  aa <- compute_g(xyz_coords,
                  xyz_coords,
                  m = m)
  ab <- compute_h(xyz_coords,
                  xyz_coords,
                  m = m)
  diag(aa) <- diag(aa) + smoothing
  aa_inv <- solve(aa)
  ag <- colSums(aa_inv)
  ag_tot <- sum(ag)
  new_sig <- as.matrix(data$signals[, data_chans])
  bb <- t(apply(new_sig,
                1,
                function(x) aa_inv %*% x))
  bb_rows <- rowSums(bb) / ag_tot
  bc <- bb - (bb_rows %*% t(ag))
  #bd <- bc %*% t(ab)
  be <- t(apply(bc,
                1,
                function(x) rowSums(x * ab)) / scaling)
  data$signals[, data_chans] <- as.data.frame(be)
  names(data$signals) <- orig_elecs
  data
}


#' Compute the h function for two sets of locations of channel locations on the
#' unit sphere.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param xyz_coords A set of electrode locations on a unit sphere.
#' @param xyz_elecs A set of electrode locations on a unit sphere.
#' @param m Interpolation constant (higher = less flexible)
#' @param lambda smoothing parameter
#' @param iter iterations for calculations
#' @importFrom pracma legendre
#' @noRd

compute_h <- function(xyz_coords,
                      xyz_elecs,
                      m = 4,
                      iter = 50) {

  EI <- 1 - sqrt((matrix(rep(xyz_coords[, 1],
                             nrow(xyz_elecs)),
                         nrow = nrow(xyz_elecs)) -
                    matrix(rep(xyz_elecs[, 1],
                               each = nrow(xyz_coords)),
                           ncol = nrow(xyz_coords))) ^ 2 +
                   (matrix(rep(xyz_coords[, 2],
                               nrow(xyz_elecs)),
                           nrow = nrow(xyz_elecs)) -
                      matrix(rep(xyz_elecs[, 2],
                                 each = nrow(xyz_coords)),
                             ncol = nrow(xyz_coords))) ^ 2 +
                   (matrix(rep(xyz_coords[, 3],
                               nrow(xyz_elecs)),
                           nrow = nrow(xyz_elecs)) -
                      matrix(rep(xyz_elecs[, 3],
                                 each = nrow(xyz_coords)),
                             ncol = nrow(xyz_coords))) ^ 2)

  dim(EI) <- c(nrow(xyz_coords), nrow(xyz_elecs))

  g <- matrix(0, ncol = ncol(EI), nrow = (EI))

  for (i in seq(1, 7)) {
    poly_xy <- pracma::legendre(i, EI)
    dim(poly_xy) <- c(i + 1,
                      nrow(EI),
                      ncol(EI))
    g <- g + ((2 * i + 1) / (i ^ m * (i + 1) ^ m)) * poly_xy[1, , ]
  }

  h <- matrix(0, ncol = ncol(EI), nrow = (EI))

  for (i in seq(1, iter)) {
    poly_xy <- pracma::legendre(i, EI)
    dim(poly_xy) <- c(i + 1,
                      nrow(EI),
                      ncol(EI))
    h <- h + ((-2 * i - 1) / (i ^ (m - 1) * (i + 1) ^ (m - 1))) * poly_xy[1, , ]
  }

  g <- g / (4 * pi)
  h <- -h / (4 * pi)
}
