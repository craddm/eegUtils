#' Create chan_info structure
#'
#' @noRd

create_chans <- function(x) {
  data.frame(chan_no)
}


#' import channel locatoins from various file formats
#' Currently only ASA .elc format is supported.
#' @param file_name name and full path of file to be loaded
#' @noRd
import_chans <- function(file_name) {
  file_type <- tools::file_ext(file_name)
  if (file_type == ".elc") {
    chan_locs <- import_elc(file_name)
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
  pos <- do.call("rbind", pos)
  final_locs <- data.frame(electrode = labs, pos3)
  names(final_locs) <- c("electrode", "x", "y", "z")
  final_locs
}


#' Convert 3D Cartesian co-ordinates to spherical
#'
#' @param x X co-ordinates
#' @param y Y co-ordinates
#' @param z z co-ordinates
#' @noRd
cart_to_sph <- function(x, y, z) {

  hypo <- sqrt(abs(x) ^2 + abs(y) ^2)
  radius <- sqrt(abs(hypo) ^2 + abs(z) ^2) # spherical radius
  phi <- atan2(z, hypo) / pi * 180 # spherical phi in degrees
  theta <- atan2(y, x) / pi * 180 # spherical theta in degrees
  data.frame(sph_radius = radius, sph_phi = phi, sph_theta = theta)
}

#' Convert 3D Cartesian co-ordinates to polar co-ordinates
#'
#' @param x X co-ordinates
#' @param y Y co-ordinates
#' @noRd

cart_to_pol <- function(x, y) {
  theta <- atan2(y, x)
  radius <- sqrt(abs(x) ^2 + abs(y) ^2)
  data.frame(pol_theta = theta, pol_radius = radius)
}
