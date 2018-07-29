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
#' @references [1] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
#'       (1989). Spherical splines for scalp potential and current
#'       density mapping. Electroencephalography and Clinical
#'     Neurophysiology, 72, 184-187
#'   [2] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
#'      (1990). Corrigenda EEG 02274. Electroencephalography and
#'      Clinical Neurophysiology, 76, 565
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

  # Note - this checks for old-style electrode locations first, then for new style
  check_ci_str(data$chan_info)

  if (all(c("cart_x", "cart_y", "cart_z") %in% names(data$chan_info))) {
    xyz_coords <- data$chan_info[, c("cart_x", "cart_y", "cart_z")]
    #normalise to unit sphere
    rads <- sqrt(rowSums(xyz_coords ^ 2))
    #rads <- sqrt(xyz_coords$cart_x ^ 2 + xyz_coords$cart_y ^ 2 + xyz_coords$cart_z ^ 2)
    xyz_coords <- xyz_coords / rads
  } else {
    xyz_coords <- sph_to_cart(data$chan_info$sph_theta / 180 * pi,
                              data$chan_info$sph_phi / 180 * pi,
                              1)
  }

  bad_select <- toupper(data$chan_info$electrode) %in% toupper(bad_elecs)
  xyz_bad <- xyz_coords[bad_select, ]
  xyz_good <- xyz_coords[!bad_select, ]

  if (nrow(xyz_bad) == 0) {
    return(data)
  }
  sigs_select <- toupper(names(data$signals)) %in%
    toupper(data$chan_info$electrode)
  bad_cols <- toupper(names(data$signals)) %in% toupper(bad_elecs)
  final_cols <- sigs_select & !bad_cols
  weights <- spheric_spline(xyz_good,
                             xyz_coords)
  new_w <- weights[bad_select, ]
  new_chans <- new_w %*% t(data$signals[, final_cols])
  data$signals[, bad_cols] <- t(new_chans)
  data
}


#' Calculate a spherical spline smooth for interpolation of electrodes
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param good_elecs Electrode positions that do not need interpolation
#' @param bad_elecs Electrode positions to be interpolated
#' @keywords internal

spheric_spline <- function(good_elecs,
                           all_elecs) {

  lec <- compute_g(good_elecs,
                   good_elecs)

  diag(lec) <- diag(lec) + 1e-05

  ph <- compute_g(all_elecs,
                  good_elecs)

  g_dims <- dim(lec)
  lec <- cbind(lec, 1)
  lec <- rbind(lec, 1)
  lec[g_dims[1] + 1, g_dims[2] + 1] <- 0
  invC <- MASS::ginv(lec, tol = 0)
  W <- cbind(ph, 1)
  W <- W %*% invC[ , 1:g_dims[1]]
  W
}


#' Calculate current source densities
#'
#'
#'
#' @param data \code{eeg_data} object
#' @param m smoothing constraint (higher = more rigid)
#' @param smoothing lambda constant
#' @param scaling Default scaling (1) is uV / m^2.
#' @keywords internal

compute_csd <- function(data,
                        m = 4,
                        smoothing = 1e-05,
                        scaling = 1) {

  orig_elecs <- names(data$signals)
  data_chans <- toupper(names(data$signals)) %in% toupper(data$chan_info$electrode)
  scaling <- scaling * scaling

  if (any(!data_chans)) {
    stop("No channel information found for ",
         names(data$signals)[!data_chans],
         ". Either remove channel or add channel info.")
  }

  # Use average reference
  data <- reref_eeg(data)

  if (all(c("cart_x", "cart_y", "cart_z") %in% names(data$chan_info))) {
    xyz_coords <- data$chan_info[, c("cart_x", "cart_y", "cart_z")]
    #normalise to unit sphere
    rads <- sqrt(xyz_coords$cart_x ^ 2 + xyz_coords$cart_y ^ 2 + xyz_coords$cart_z ^ 2)
    xyz_coords <- xyz_coords / rads
  } else {
    xyz_coords <- sph_to_cart(data$chan_info$sph_theta / 180 * pi,
                              data$chan_info$sph_phi / 180 * pi,
                              1)
  }

  g_mat <- compute_g(xyz_coords,
                  xyz_coords,
                  m = m,
                  iter = 50)
  h_mat <- compute_h(xyz_coords,
                  xyz_coords,
                  m = m,
                  iter = 50)

  diag(g_mat) <- diag(g_mat) + smoothing
  g_inv <- solve(g_mat)
  ag <- colSums(g_inv)
  ag_tot <- sum(ag)
  new_sig <- as.matrix(data$signals[, data_chans])
  bb <- t(apply(new_sig,
                1,
                function(x) g_inv %*% x))
  bb_rows <- rowSums(bb) / ag_tot
  bc <- bb - (bb_rows %*% t(ag))
  be <- t(apply(bc,
                1,
                function(x) colSums(x * h_mat)) / scaling)
  data$signals[, data_chans] <- as.data.frame(be)
  names(data$signals) <- orig_elecs
  data
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
#' @keywords internal

compute_g <- function(xyz_coords,
                      xyz_elecs,
                      m = 4,
                      iter = 7) {

  xyz_coords <- as.matrix(xyz_coords)
  xyz_elecs <- as.matrix(xyz_elecs)

  EI <- matrix(NA,
               nrow = nrow(xyz_coords),
               ncol = nrow(xyz_coords))
  EI <- xyz_coords %*% t(xyz_elecs) # cosine similarity

  g <- matrix(0, ncol = ncol(EI), nrow = nrow(EI))

  for (i in seq(1, iter)) {
    suppressWarnings(
      poly_xy <- pracma::legendre(i, EI)
    )
    dim(poly_xy) <- c(i + 1,
                      nrow(EI),
                      ncol(EI))
    g <- g + ((2 * i + 1) / (i ^ m * (i + 1) ^ m)) * poly_xy[1, , ]
  }
  g <- g / 4 / pi
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
#' @keywords internal

compute_h <- function(xyz_coords,
                      xyz_elecs,
                      m = 4,
                      iter = 50) {

  xyz_coords <- as.matrix(xyz_coords)
  xyz_elecs <- as.matrix(xyz_elecs)
  EI <- xyz_coords %*% t(xyz_elecs)

  h <- matrix(0,
              ncol = ncol(EI),
              nrow = nrow(EI))

  for (i in seq(1, iter)) {
    suppressWarnings(
      poly_xy <- pracma::legendre(i, EI)
    )
    dim(poly_xy) <- c(i + 1,
                      nrow(EI),
                      ncol(EI))
    h <- h + ((-2 * i - 1) / (i ^ (m - 1) * (i + 1) ^ (m - 1))) * poly_xy[1, , ]
  }

  h <- -h / 4 / pi
}
