#' Channel interpolation
#'
#' Interpolate EEG channels using a spherical spline (Perrin et al., 1989; 1990). The
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
#' [2] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
#'      (1990). Corrigenda EEG 02274. Electroencephalography and
#'      Clinical Neurophysiology, 76, 565
#' @export

interp_elecs <- function(data, bad_elecs, ...) {
  UseMethod("interp_elecs", data)
}

#' @export
interp_elecs.default <- function(data, bad_elecs, ...) {
  stop("Not implemented for objects of class", class(data))
}

#' @describeIn interp_elecs Interpolate EEG channel(s)
#' @export
interp_elecs.eeg_data <- function(data,
                                  bad_elecs,
                                  ...) {

  if (is.null(data$chan_info)) {
    stop("No channel locations found.")
  }

  bads_check <- bad_elecs %in% data$chan_info$electrode

  if (!all(bads_check)) {
    stop("No specified channels found.")
  }

  if (any(!bads_check)) {
    warning("Electrode(s) not found: ",
            paste0(bad_elecs[!bads_check],
                   collapse = " "))
    bad_elecs <- bad_elecs[bads_check]
  }

  data$chan_info <- validate_channels(data$chan_info,
                                      channel_names(data))
  missing_coords <- apply(is.na(data$chan_info), 1, any)
  missing_chans <- names(data$signals)[missing_coords]

  if (any(missing_coords)) {
    message("Coords missing for electrodes ",
            paste0(missing_chans,
                   collapse = " "))
  }

  chan_info <- data$chan_info[!missing_coords, ]

  # Note - this checks for old-style electrode locations first, then for new
  # style
  check_ci_str(chan_info)

  if (all(c("cart_x", "cart_y", "cart_z") %in% names(chan_info))) {
    xyz_coords <- chan_info[, c("cart_x", "cart_y", "cart_z")]

    #normalise to unit sphere
    #rads <- sqrt(rowSums(xyz_coords ^ 2))
    #xyz_coords <- xyz_coords / rads
    xyz_coords <- norm_sphere(xyz_coords)
  } else {
    xyz_coords <- sph_to_cart(chan_info$theta,
                              chan_info$phi,
                              1)
  }

  bad_select <- toupper(chan_info$electrode) %in% toupper(bad_elecs)
  xyz_bad <- xyz_coords[bad_select, ]
  xyz_good <- xyz_coords[!bad_select, ]

  if (nrow(xyz_bad) == 0) {
    return(data)
  }

  sigs_select <- toupper(names(data$signals)) %in%
    toupper(chan_info$electrode)

  bad_cols <- toupper(names(data$signals)) %in% toupper(bad_elecs)
  final_cols <- sigs_select & !bad_cols

  weights <- spheric_spline(xyz_good,
                             xyz_coords)

  data$signals <- interp_chans(data$signals,
                          bad_elecs,
                          missing_coords = missing_coords,
                          weights)
  data
}


#' Calculate a spherical spline smooth for interpolation of electrodes
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param good_elecs Electrode positions that do not need interpolation
#' @param all_elecs Electrode positions including those to be interpolated
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

#' Interpolate channels
#' @param .data Channel data containing all data
#' @param bad_chans Vector of names of bad channels
#' @param missing_coords Logical vector indicating any channels in the data that
#'   had no associated coordinates
#' @param weights Spherical spline weights for interpolation
#' @keywords internal
interp_chans <- function(.data,
                         bad_chans,
                         missing_coords = FALSE,
                         weights) {

  bad_cols <- (toupper(names(.data)) %in% toupper(bad_chans)) | missing_coords
  weight_rows <- names(.data[, !missing_coords]) %in% bad_chans
  new_chans <- weights[weight_rows, , drop = TRUE] %*% t(.data[ , !bad_cols])
  .data[, bad_chans] <- t(new_chans)
  .data
}

#' Convert to Current Source Density
#'
#' Convert an \code{eeg_data} or \code{eeg_epochs} object to using Current
#' Source Densities. This command uses a spherical spline algorithm (Perrin et
#' al., 1989) to compute scalp surface Laplacian/current source density
#' estimates of scalp potentials, a reference-free measure of electrical
#' activity that emphasises more local spatial features
#'
#' @param data \code{eeg_data} or \code{eeg_epochs} object
#' @param ... Other parameters
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @references [1] Perrin, F., Pernier, J., Bertrand, O., Echallier, J.F.
#'   (1989). Spherical splines for scalp potential and current density mapping.
#'   Electroencephalography and Clinical Neurophysiology, 72(2), 184-187. PMID:
#'   2464490 [2] Kayser, J., Tenke, C.E. (2006). Principal components analysis
#'   of Laplacian waveforms as a generic method for identifying ERP generator
#'   patterns: I. Evaluation with auditory oddball tasks. Clinical
#'   Neurophysiology, 117(2), 348-368. [3] Kayser, J., Tenke, C.E. (2015).
#'   Issues and considerations for using the scalp surface Laplacian in EEG/ERP
#'   research: A tutorial review. International Journal of Psycholphysiology,
#'   97(3), 189-209
#' @examples
#' csd_epochs <- compute_csd(demo_epochs)
#' plot_butterfly(csd_epochs)
#' @export
compute_csd <- function(data,
                        ...) {
  UseMethod("compute_csd", data)
}

#' @describeIn compute_csd Default method to detect unknown classes.
#' @export
compute_csd.default <- function(data,
                                ...) {
  warning("compute_csd not implemented for objects of class", class(data))
}

#' @param m smoothing constraint (higher = more rigid splines)
#' @param smoothing lambda constant. Added to the Defaults to 1e-05
#' @param scaling Default scaling (1) is uV / m^2. Note that this depends on the
#'   units of the electrode co-ordinates.
#' @describeIn compute_csd Transform eeg_data to CSD
#' @export
compute_csd.eeg_data <- function(data,
                                 m = 4,
                                 smoothing = 1e-05,
                                 scaling = 1,
                                 ...) {
  convert_to_csd(data, m, smoothing, scaling)
}

#' @describeIn compute_csd Transform eeg_data to CSD
#' @export
compute_csd.eeg_epochs <- function(data,
                                 m = 4,
                                 smoothing = 1e-05,
                                 scaling = 1,
                                 ...) {
  convert_to_csd(data, m, smoothing, scaling)
}

#' Calculate current source densities
#'
#' @param data \code{eeg_data} object
#' @param m smoothing constraint (higher = more rigid)
#' @param smoothing lambda constant
#' @param scaling Default scaling (1) is uV / m^2.
#' @keywords internal

convert_to_csd <- function(data,
                           m = 4,
                           smoothing = 1e-05,
                           scaling = 1) {

  orig_elecs <- names(data$signals)
  data_chans <- toupper(names(data$signals)) %in% toupper(data$chan_info$electrode)
  scaling <- scaling * scaling

  if (any(!data_chans)) {
    stop("No channel information found for ",
         paste0(names(data$signals)[!data_chans],
                collapse = " "),
         ". Either remove channel or add channel info.")
  }

  missing_coords <- apply(data$chan_info, 1, function(x) any(is.na(x)))

  if (any(missing_coords)) {
    stop("No coordinates for ",
         paste0(data$chan_info$electrode[missing_coords],
                collapse = " ")
         )
  }

  # Convert data to average reference
  data <- reref_eeg(data)

  if (all(c("cart_x", "cart_y", "cart_z") %in% names(data$chan_info))) {
    xyz_coords <- data$chan_info[, c("cart_x", "cart_y", "cart_z")]
    #normalise to unit sphere
    #rads <- sqrt(xyz_coords$cart_x ^ 2 + xyz_coords$cart_y ^ 2 + xyz_coords$cart_z ^ 2)
    #xyz_coords <- xyz_coords / rads
    xyz_coords <- norm_sphere(xyz_coords)
  } else {
    xyz_coords <- sph_to_cart(data$chan_info$theta,
                              data$chan_info$phi,
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
  sums_g <- colSums(g_inv)
  total_g <- sum(sums_g)
  new_sig <- as.matrix(data$signals[, data_chans])
  bb <- t(apply(new_sig,
                1,
                function(x) g_inv %*% x))
  bb_rows <- rowSums(bb) / total_g
  bc <- bb - (bb_rows %*% t(sums_g))
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

  gg <- 1:iter
  gg <- (2 * gg + 1) / (gg ^ m *(gg + 1) ^ m)
  legpoly <- matrix(0,
                    nrow = length(c(EI)),
                    ncol = iter)

  for (i in seq(1, iter)) {
    suppressWarnings(
      poly_xy <- pracma::legendre(i, EI)
    )
    legpoly[, i] <- t(poly_xy[1, ])
  }

  g <- sweep(legpoly, 2, gg, "*")#g + gg
  g <- rowSums(g)
  g <- g / 4 / pi
  dim(g) <- c(nrow(EI), ncol(EI))
  g
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

  # vectorize this too, just like compute_g!
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
  h
}
