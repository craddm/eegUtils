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
#' @noRd

import_elc <- function(file_name) {
  raw_locs <- readLines(file_name, n = -1)
  n_elecs <- grep("NumberPositions", raw_locs)
  n_elecs <- as.numeric(unlist(strsplit(raw_locs[n_elecs], "\t"))[2])
  pos_loc <- grep("^Positions", raw_locs)
  pos <- raw_locs[seq(pos_loc + 1, pos_loc + n_elecs)]
  labs_loc <- grep("Labels", raw_locs)
  labs <- raw_locs[seq(labs_loc + 1, labs_loc + n_elecs)]

  pos <- strsplit(pos, " ")
  pos <- lapply(pos, function(x) as.numeric(x[!x == ""]))
  pos <- as.data.frame(do.call("rbind", pos))
  sph_pos <- cart_to_sph(pos[, 1], pos[, 2], pos[, 3])
  topo_pos <- sph_to_topo(sph_pos[, 2], sph_pos[, 3])

  pol_pos <- cart_to_pol(pos[, 1], pos[, 2])
  names(pos) <- c("cart_x", "cart_y", "cart_z")
  final_locs <- data.frame(electrode = labs,
                           pos,
                           sph_pos,
                           pol_pos,
                           topo_pos)
  final_locs$x <- final_locs$radius * cos(final_locs$angle / 180 * pi)
  final_locs$y <- final_locs$radius * sin(final_locs$angle / 180 * pi)
  final_locs
}

#' Convert 3D Cartesian co-ordinates to spherical
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param x X co-ordinates
#' @param y Y co-ordinates
#' @param z z co-ordinates
#' @return Data frame with entries "sph_radius", "sph_phi" (in degrees), "sph_theta" (in degrees).
#' @noRd

cart_to_sph <- function(x, y, z) {

  hypo <- sqrt(abs(x) ^ 2 + abs(y) ^ 2)
  radius <- sqrt(abs(hypo) ^ 2 + abs(z) ^ 2) # spherical radius
  phi <- atan2(z, hypo) / pi * 180 # spherical phi in degrees
  theta <- atan2(y, x) / pi * 180 # spherical theta in degrees
  data.frame(sph_radius = radius, sph_phi = phi, sph_theta = theta)
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

#' Convert polar to spherical coordinates
#'
#' Hard-coded to a radius of 85 mm (as in BESA).
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param theta Azimuth from polar co-ordinates (theta in supplied electrode
#'   locations)
#' @param phi Elevation from polar co-ordinates (radius in supplied electrode
#'   locations)
#' @noRd

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
#' @param phi Phi
#' @param theta Theta
#' @noRd
sph_to_topo <- function(phi, theta) {
  angle <- -theta
  radius <- 0.5 - phi / 180
  data.frame(angle, radius)
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
#' @import dplyr
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
              paste(dplyr::pull(unique(data[, electrode]))[!elecs], sep = ","))
    } else if (!any(elecs)) {
      stop("No matching electrodes found.")
    }
  } else {
    elecs <-
      unique(data[, electrode]) %in% dplyr::pull(electrodeLocs[, electrode])
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
#' @import dplyr
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
  missing_els <- elec_names %in% electrodeLocs$electrode

  if (!any(matched_els)) {
    stop("No matching electrodes found.")
  } else if (!all(matched_els)) {
    message(cat("Electrodes not found:", names(data$signals)[!missing_els]))
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

plot_electrodes.default <- function(data, interact = FALSE) {

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
#' @noRd

montage_check <- function(montage) {
  if (identical(montage, "biosemi64alpha")) {
    elocs <- orig_locs
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
