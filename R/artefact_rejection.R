#' FASTER EEG artefact rejection
#'
#' Whelan et al (2011). Not yet implemented.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_epochs}
#' @param ... Parameters passed to FASTER

eeg_FASTER <- function(data, ...) {
  UseMethod("eeg_FASTER", data)
}

#' @describeIn eeg_FASTER Run FASTER on \code{eeg_epochs}
eeg_FASTER.eeg_epochs <- function(data, ...) {


}

#' Simple absolute value thresholding
#'
#' Reject data based on a simple absolute threshold. This marks any
#' timepoint from any electrode.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class \code{eeg_data} or \code{eeg_epochs}.
#' @param threshold In microvolts. If one value is supplied, it will be treated
#'   as a +- value.
#' @param reject If TRUE, remove marked data immediately, otherwise mark for
#'   inspection/rejection. Defaults to FALSE.
#' @param ... Other arguments passed to eeg_ar_thresh
#' @export

eeg_ar_thresh <- function(data, threshold, reject = FALSE, ...) {
  UseMethod("eeg_ar_thresh", data)
}

#' @describeIn eeg_ar_thresh Reject data using a simple threshold.
eeg_ar_thresh.eeg_data <- function(data, threshold, reject = FALSE, ...) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- data$signals > max(threshold) |
    data$signals < min(threshold)

  if (reject) {
    crossed_thresh <- rowSums(crossed_thresh) == 0
    data$timings <- data$timings[crossed_thresh, ]
    data$signals <- data$signals[crossed_thresh, ]
    data$events <- data$events[data$events$event_time %in% data$timings$time, ]
    data$reference$ref_data <- data$reference$ref_data[crossed_thresh, ]
  } else {
    data$reject <- crossed_thresh
  }
  data
}

#' @describeIn eeg_ar_thresh Reject data using a simple threshold.
eeg_ar_thresh.eeg_epochs <- function(data, threshold, reject = FALSE, ...) {

  if (length(threshold) == 1) {
    threshold <- c(threshold, -threshold)
  }

  crossed_thresh <- data$signals > max(threshold) | data$signals < min(threshold)

  if (reject) {
    crossed_thresh <- rowSums(crossed_thresh) == 1
    rej_epochs <- unique(data$timings$epoch[crossed_thresh])
    data <- select_epochs(data, rej_epochs, keep = FALSE) # consider creating select_timerange vs select_timepoints
  } else {
    data$reject <- rej_epochs
  }
  data
}

#' Channel statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param ... Other parameters passed to the functions.


channel_stats <- function(data, ...) {
  UseMethod("channel_stats", data)
}

#' @describeIn channel_stats Calculate channel statistics for \code{eeg_data} objects.
channel_stats.eeg_data <- function(data, ...) {
  chan_means <- colMeans(data$signals)
  chan_sds <- apply(data$signals, 2, sd)
  chan_var <- apply(data$signals, 2, var)
  chan_kurt <- apply(data$signals, 2, kurtosis)
}

#' Epoch statistics
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data as a \code{eeg_data} or \code{eeg_epochs} object.
#' @param ... Other parameters passed to the functions.

epoch_stats <- function(data, ...) {
  UseMethod("epoch_stats", data)
}

#' @describeIn epoch_stats Calculate statistics for each epoch.
epoch_stats.eeg_epochs <- function(data, ...) {
  epoch_nos <- unique(data$timings$epoch)
  data <- split(data$signals, data$timings$epoch)
  epoch_means <- lapply(data, colMeans)
  epoch_means <- as.data.frame(do.call("rbind", epoch_means))
  epoch_means$epoch_nos <- epoch_nos
  epoch_means
}

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
interp_elecs.eeg_data <- function(data, bad_elecs, ...) {

  if (is.null(data$chan_info)) {
    stop("No channel locations found.")
  }
  xyz_coords <- pol_to_sph(data$chan_info$theta, data$chan_info$radius)
  rads <- sqrt(xyz_coords$x ^ 2 + xyz_coords$y ^ 2 + xyz_coords$z ^ 2)
  xyz_coords <- xyz_coords / rads

  bad_select <- toupper(data$chan_info$electrode) %in% toupper(bad_elecs)
  xyz_bad <- xyz_coords[bad_select, ]
  xyz_good <- xyz_coords[!bad_select, ]

  sigs_select <- toupper(names(data$signals)) %in%
    toupper(data$chan_info$electrode)

  sigs <- data$signals[sigs_select]

  bad_cols <- toupper(names(sigs)) %in% toupper(bad_elecs)

  mm <- spheric_spline(xyz_good, xyz_bad, sigs[!bad_cols])

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

  C <- as.matrix(data) %*% t(MASS::ginv(rbind(lec, 1)))

  #C <- C[-ncol(C)]
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
#' @importFrom pracma legendre
#' @noRd

compute_g <- function(xyz_coords, xyz_elecs) {

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
    dim(poly_xy) <- c(i + 1, nrow(EI), ncol(EI))
    g <- g + ((2 * i + 1) / (i ^ 4 * (i + 1) ^ 4)) * poly_xy[1, , ]
  }

  g <- g / (4 * pi)

}


#' Calculate kurtosis
#'
#' @param data Data to calculate kurtosis for
#' @noRd

kurtosis <- function(data) {
  m4 <- mean((data - mean(data)) ^ 4)
  kurt <- m4 / (sd(data) ^ 4) - 3
  kurt
}
