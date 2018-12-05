#' Import channel locations from various file formats
#'
#' Currently only ASA .elc format with Cartesian x-y-z coordinates is supported.
#'
#' @param file_name name and full path of file to be loaded
#' @export

import_chans <- function(file_name) {
  file_type <- tools::file_ext(file_name)
  if (file_type == "elc") {
    chan_locs <- import_elc(file_name)
  } else {
    stop("File type ", file_type, " is unknown.")
  }
}

#' Import ASA .elc electrode location files
#'
#' Loads and process ASA electrode locations.
#'
#' @param file_name file name
#' @keywords internal

import_elc <- function(file_name) {
  raw_locs <- readLines(file_name,
                        n = -1)
  n_elecs <- grep("NumberPositions",
                  raw_locs)
  n_elecs <- as.numeric(unlist(strsplit(raw_locs[n_elecs], "\t"))[2])
  pos_loc <- grep("^Positions", raw_locs)
  pos <- raw_locs[seq(pos_loc + 1,
                      pos_loc + n_elecs)]
  labs_loc <- grep("Labels", raw_locs)
  labs <- raw_locs[seq(labs_loc + 1,
                       labs_loc + n_elecs)]

  pos <- strsplit(pos, " ")
  pos <- lapply(pos,
                function(x) as.numeric(x[!x == ""]))
  pos <- as.data.frame(do.call("rbind", pos))
  sph_pos <- cart_to_sph(pos[, 1],
                         pos[, 2],
                         pos[, 3])

  topo_pos <- sph_to_topo(phi = sph_pos[, 2],
                          theta = sph_pos[, 3])
  pol_coords <- cart_to_pol(pos[, 1],
                            pos[, 2])
  names(pos) <- c("cart_x",
                  "cart_y",
                  "cart_z")
  final_locs <- data.frame(electrode = labs,
                           pos,
                           sph_pos,
                           pol_coords,
                           topo_pos)#,
  final_locs <- cbind(final_locs,
                      topo_norm(final_locs$angle,
                                final_locs$radius))

  #final_locs$x <- final_locs$radius * cos(final_locs$angle / 180 * pi)
  #final_locs$y <- final_locs$radius * sin(-final_locs$angle / 180 * pi)
  tibble::as.tibble(final_locs)
}

#' Convert 3D Cartesian co-ordinates to spherical
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x X co-ordinates
#' @param y Y co-ordinates
#' @param z Z co-ordinates
#' @return Data frame with entries "sph_radius",
#'  "sph_phi" (in degrees),
#'   "sph_theta" (in degrees).
#' @keywords internal

cart_to_sph <- function(x, y, z) {

  hypo <- sqrt(abs(x) ^ 2 + abs(y) ^ 2)
  radius <- sqrt(abs(hypo) ^ 2 + abs(z) ^ 2) # spherical radius
  phi <- atan2(z, hypo)  * 180 / pi# spherical phi in degrees
  theta <- atan2(y, x)  * 180 / pi# spherical theta in degrees
  data.frame(sph_radius = radius,
             sph_phi = phi,
             sph_theta = theta)
}

#' Convert 3D Cartesian co-ordinates to polar co-ordinates
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x X co-ordinates
#' @param y Y co-ordinates
#' @noRd

cart_to_pol <- function(x, y) {
  theta <- atan2(y, x) / pi * 180
  radius <- sqrt(abs(x) ^ 2 + abs(y) ^ 2)
  data.frame(pol_theta = theta, pol_radius = radius)
}

#' Convert EEGLAB polar to spherical coordinates
#'
#' Hard-coded to a radius of 85 mm (as in BESA).
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param theta Azimuth from polar co-ordinates (theta in supplied electrode
#'   locations)
#' @param phi Elevation from polar co-ordinates (radius in supplied electrode
#'   locations)
#' @keywords internal

pol_to_sph <- function(theta, phi) {

  theta <- -theta
  phi <- (0.5 - phi) * 180
  sph_theta <- theta / 180 * pi
  sph_phi <- phi / 180 * pi
  radius <- 85

  z <- radius * sin(sph_phi)
  z_cos <- radius * cos(sph_phi)
  x <- z_cos * cos(sph_theta)
  y <- z_cos * sin(sph_theta)
  data.frame(x, y, z)
}

#' Convert spherical to topographical co-ordinates
#'
#' Expects input in degrees
#'
#' @param phi Phi
#' @param theta Theta
#' @keywords internal
sph_to_topo <- function(theta, phi) {

  angle <- -theta
  radius <- 0.5 - phi / 180
  data.frame(angle, radius)
}

#' Convert spherical to cartesian 3d
#'
#' Note that phi follows EEGLAB conventions (i.e. needs correcting to 90 - phi in some cases)
#' @param theta should be in degrees
#' @param phi should be in degrees
#' @param r should be in degrees
#' @keywords internal
sph_to_cart <- function(theta, phi, radius) {
  z <- radius * sin(phi * pi / 180)
  x <- radius * cos(phi * pi / 180) * cos(theta * pi / 180)
  y <- radius * cos(phi * pi / 180) * sin(theta * pi / 180)
  data.frame(cart_x = x, cart_y = y, cart_z = z)
}

#' Convert topographical 2d to cartesian 2d
#'
#' Expects input to be in degrees
#'
#' @param angle Angle
#' @param radius Radius
#' @keywords internal

topo_norm <- function(angle, radius) {
  x <- radius * cos(angle * pi / 180)
  y <- radius * sin(-angle * pi / 180)
  data.frame(x, y)
}

#' Rotate channel locations
#'
#' @param chan_info channel information structure
#' @param degrees degrees by which to rotate
#' @keywords internal

rotate_angle <- function(chan_info, degrees) {

  degrees <- degrees * pi / 180
  if ("CZ" %in% chan_info$electrode) {
    cent_x <- chan_info[toupper(chan_info$electrode) == "Cz", ]$x
    cent_y <- chan_info[toupper(chan_info$electrode) == "Cz", ]$y
  } else {
    cent_x <- 0
    cent_y <- 0
  }

  chan_info$x <- chan_info$x - cent_x
  chan_info$y <- chan_info$y - cent_y
  rot_x <- cent_x + cos(degrees) * chan_info$x - sin(degrees) * chan_info$y
  rot_y <- cent_y + sin(degrees) * chan_info$x + cos(degrees) * chan_info$y
  chan_info$x <- rot_x
  chan_info$y <- rot_y
  chan_info
}

#' Rotate spherical coordinates and recalculate others
#'
#' @param chan_info channel information structure
#' @param degrees degrees by which to rotate elecs
#' @keywords internal
rotate_sph <- function(chan_info, degrees) {
  chan_info$sph_theta <- chan_info$sph_theta + degrees
  chan_info$sph_theta <- ifelse(chan_info$sph_theta > 180,
                            chan_info$sph_theta - 360,
                            chan_info$sph_theta)
  chan_info$sph_theta <- ifelse(chan_info$sph_theta < -180,
                                chan_info$sph_theta + 360,
                                chan_info$sph_theta)
  topo_pos <- sph_to_topo(theta = chan_info$sph_theta,
                          phi = chan_info$sph_phi)
  chan_info$angle <- topo_pos[, 1]
  chan_info$radius <- topo_pos[, 2]
  cart_sph <- pol_to_sph(chan_info$angle,
                         phi = chan_info$radius)
  chan_info$cart_x <- cart_sph[, 1]
  chan_info$cart_y <- cart_sph[, 2]
  chan_info$cart_z <- cart_sph[, 3]
  chan_info$x <- chan_info$radius * cos(chan_info$angle / 180 * pi)
  chan_info$y <- chan_info$radius * sin(chan_info$angle / 180 * pi)
  chan_info
}


#' Flip x-axis coords
#'
#' @param chan_info chan-info structure
#' @keywords internal

flip_x <- function(chan_info) {
  chan_info$cart_x <- chan_info$cart_x * -1
  chan_info$cart_y <- chan_info$cart_y * -1
  chan_info$x <- chan_info$x * -1
  chan_info$angle <- chan_info$angle * -1
  chan_info$sph_theta <- chan_info$sph_theta * -1
  chan_info
}
#' Get standard electrode locations
#'
#' Joins standard electrode locations to EEG data from eegUtils internal data.
#'
#' @param data An EEG dataset.
#' @param ... Parameters passed to S3 methods.
#' @export

electrode_locations <- function(data, ...) {
  UseMethod("electrode_locations")
}

#' @param electrode The column name containing electrode names in data.
#'   (Defaults to "electrode").
#' @param drop Should electrodes in \code{data} for which default locations are
#'   not available be dropped? (Defaults to FALSE).
#' @param plot Plot obtained electrode locations.
#' @param montage Name of an existing montage set. Defaults to NULL; (currently
#'   only 'biosemi64alpha' available other than default 10/20 system)
#' @importFrom dplyr inner_join pull left_join distinct
#' @import ggplot2
#' @importFrom tibble is.tibble
#' @describeIn electrode_locations Adds standard locations to a data frame in
#'   long format
#' @return A tibble (or data.frame), or ggplot2 object if \code{plot = TRUE}.
#' @export

electrode_locations.data.frame <- function(data,
                                           electrode = "electrode",
                                           drop = FALSE,
                                           plot = FALSE,
                                           montage = NULL, ...) {

  #if a montage supplied, check if it matches known montages
  if (!is.null(montage)) {
    electrodeLocs <- montage_check(montage)
  }

  data[, electrode] <- toupper(data[[electrode]])
  electrodeLocs[, electrode] <- toupper(electrodeLocs[[electrode]])

  if (tibble::is.tibble(data)) {
    elecs <-
      dplyr::pull(unique(data[, electrode])) %in%
      dplyr::pull(electrodeLocs[, electrode])

    if (!all(elecs)) {
      message("Electrodes not found: ",
              paste(unique(data[, electrode])[!elecs],
                    sep = ","))
    } else if (!any(elecs)) {
      stop("No matching electrodes found.")
    }
  } else {
    elecs <-
      unique(data[, electrode]) %in% electrodeLocs[, electrode,
                                                   drop = TRUE]
    if (!all(elecs)) {
      message("Electrodes not found: ",
              paste(unique(data[, electrode])[!elecs], sep = ","))
    } else if (!any(elecs)) {
      stop("No matching electrodes found.")
    }

  }

  if (drop) {
    data <- dplyr::inner_join(data, electrodeLocs, by = electrode)
  } else {
    data <- dplyr::left_join(data, electrodeLocs, by = electrode)
  }

  if (plot) {
    plotdata <- dplyr::distinct(data, x, y, electrode)
    p <- ggplot2::ggplot(plotdata, aes(x, y)) +
      geom_label(aes(label = electrode))
    return(p)
  } else {
    return(data)
  }
}

#' @param overwrite Overwrite existing channel info. Defaults to FALSE.
#' @import ggplot2
#' @describeIn electrode_locations Adds standard locations to the chan_info field of an eeg_data object.
#' @export

electrode_locations.eeg_data <- function(data,
                                         drop = FALSE,
                                         plot = FALSE,
                                         montage = NULL,
                                         overwrite = FALSE, ...) {

  if (!is.null(data$chan_info) & !overwrite & !plot) {
    stop("Channel info already present, set overwrite to TRUE to replace.")
  }

  if (!is.null(montage)) {
    electrodeLocs <- montage_check(montage)
  }

  elec_names <- toupper(names(data$signals))
  electrodeLocs$electrode <- toupper(electrodeLocs$electrode)

  matched_els <- electrodeLocs$electrode %in% elec_names
  missing_els <- !elec_names %in% electrodeLocs$electrode

  if (!any(matched_els)) {
    stop("No matching electrodes found.")
  } else if (any(missing_els)) {
    message(paste("Electrodes not found:", names(data$signals)[missing_els]))
  }

  data$chan_info <- electrodeLocs[matched_els, ]

  if (drop) {
    data$signals[matched_els]
  }

  if (plot) {
    p <- ggplot2::ggplot(data$chan_info, aes(x, y)) +
      geom_label(aes(label = electrode))
    return(p)
  } else {
    return(data)
  }
}

#' Plot electrode locations
#'
#' Produces either a 2D plot of the electrode locations or an interactive plot
#' of electrode locations in 3D space.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data Data with associated electrode locations to be plotted.
#' @param interact Choose 2D cartesian layout, or, if set to TRUE, an
#'   interactive 3D plot of electrode locations. Defaults to FALSE.
#' @export

plot_electrodes <- function(data, interact = FALSE) {
  UseMethod("plot_electrodes", data)
}

#' @import ggplot2
#' @describeIn plot_electrodes generic plot electrodes function
#' @export

plot_electrodes.default <- function(data,
                                    interact = FALSE) {

  if ("electrode" %in% names(data)) {
    data <- data.frame(electrode = unique(data$electrode))
    data <- electrode_locations(data)

    if (interact) {
      if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Package \"plotly\" needed for interactive electrode plots. Please install it.",
             call. = FALSE)
      }
      plotly::plot_ly(data,
                      x = ~cart_x,
                      y = ~cart_y,
                      z = ~cart_z,
                      text = ~electrode,
                      type = "scatter3d",
                      mode = "text+markers")
    } else {
      ggplot2::ggplot(data,
                      aes(x = x,
                          y = y,
                          label = electrode)) +
        geom_text() +
        theme_minimal() +
        coord_equal()
    }
  } else {
    stop("No electrodes found.")
  }
}

#' @describeIn plot_electrodes Plot electrodes associated with an \code{eeg_data} object.
#' @export
plot_electrodes.eeg_data <- function(data,
                                     interact = FALSE) {

  if (is.null(data$chan_info)) {
    warning("Adding standard locations...")
    data <- electrode_locations(data)
  }

  if (interact) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package \"plotly\" needed for interactive electrode plots. Please install it.",
           call. = FALSE)
    }

    plotly::plot_ly(data$chan_info,
                    x = ~cart_x,
                    y = ~cart_y,
                    z = ~cart_z,
                    text = ~electrode,
                    type = "scatter3d",
                    mode = "text+markers")
  } else {
    ggplot2::ggplot(data$chan_info,
                    aes(x = x,
                        y = y,
                        label = electrode)) +
      geom_text() +
      theme_minimal() +
      coord_equal()
  }
}

#' Montage check
#'
#' @param montage Name of montage
#' @keywords internal

montage_check <- function(montage) {
  if (identical(montage, "biosemi64alpha")) {
    elocs <- merge(orig_locs["electrode"][1:64, ],
                   electrodeLocs,
                   sort = FALSE) #hacky way to translate elec names
    elocs[1:64, "electrode"] <- c(paste0("A", 1:32),
                                          paste0("B", 1:32))
  } else {
    stop("Unknown montage. Current only biosemi64alpha is available.")
  }
  elocs
}


#' Create chan_info structure
#'
#' @param chans Channel numbers
#' @param elecs Electrode names
#' @noRd

create_chans <- function(chans, elecs) {
  stopifnot(is.numeric(chans),
            is.character(elecs))
  data.frame(chan_no = chans,
             electrode = elecs)

}

empty_chans <- function() {
  data.frame(electrode = character(),
             cart_x = numeric(),
             cart_y = numeric(),
             cart_z = numeric(),
             sph_radius = numeric(),
             sph_ph = numeric(),
             sph_theta = numeric(),
             pol_theta = numeric(),
             pol_radius = numeric(),
             angle = numeric(),
             radius = numeric(),
             x = numeric(),
             y = numeric())
}
