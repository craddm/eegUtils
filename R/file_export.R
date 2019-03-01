# Export continuous data in Brain Vision Analyzer format
#' @noRd
export_bva <- function(.data,
                       filename,
                       orientation) {
  UseMethod("export_bva", .data)
}

export_bva.default <- function(.data,
                               filename,
                               orientation) {
  stop("export_bva() can currently only export continuous eeg_data objects.")
}

export_bva.eeg_epochs <- function(.data,
                                  filename,
                                  orientation = "VECTORIZED") {

  if (!(orientation %in% c("VECTORIZED", "MULTIPLEXED"))) {
    stop("Orientation must be VECTORIZED or MULTIPLEXED.")
  }

  filename <- tools::file_path_sans_ext(filename)

  write_vhdr(.data, filename, orientation)
  write_dat(.data, filename, orientation)
  write_vmrk(.data, filename)
}


write_vhdr <- function(.data,
                       filename,
                       orientation) {
  vhdr_file <- paste0(filename, ".vhdr")
  con <- file(vhdr_file, open = "w", encoding = "UTF-8")

  on.exit(close(con))
  new_header <- list()
  new_header[["Common Infos"]] <-
    list(Codepage = "UTF-8",
         DataFile = paste0(filename, ".dat"),
         MarkerFile = paste0(filename, ".vmrk"),
         DataFormat = "BINARY",
         DataOrientation=orientation,
         DataType = "TIMEDOMAIN",
         NumberOfChannels = ncol(.data$signals),
         DataPoints = nrow(.data$signals) * ncol(.data$signals),
         SamplingInterval = 1e6 / .data$srate)
  #new_header[["User Infos"]] <- ""
  new_header[["Binary Infos"]] <- list(BinaryFormat = "IEEE_FLOAT_32")
  new_header[["Channel Infos"]] <- as.list(paste(names(.data$signals),
                                                 .data$reference$ref_chans,
                                                 "",
                                                 paste0(intToUtf8(0x03BC), "V"),
                                                 sep = ","))
  names(new_header[["Channel Infos"]]) <- paste0("Ch", 1:ncol(.data$signals))
  new_header[["Coordinates"]] <- as.list(paste(.data$chan_info$sph_radius,
                                               .data$chan_info$sph_theta,
                                               .data$chan_info$sph_phi,
                                               sep = ","))
  names(new_header[["Coordinates"]]) <- paste0("Ch", 1:ncol(.data$signals))
  writeLines("Brain Vision Data Exchange Header File Version 2.0", con)
  writeLines(paste0("; Created using eegUtils http://craddm.github.io/eegUtils/"),
             con)
  writeLines("", con)

  for (i in names(new_header)) {
    writeLines(paste0("[", i,"]"), con)
    entries <-
      lapply(seq_along(new_header[[i]]),
             function(x) paste0(names(new_header[[i]][x]),
                                "=",
                                new_header[[i]][x]))
    lapply(entries,
           writeLines,
           con)
    writeLines("", con)
  }
  writeLines("", con)
}

write_dat <- function(.data,
                      filename,
                      orientation) {
  vdat_file <- paste0(filename, ".dat")
  con <- file(vdat_file, open = "wb")
  on.exit(close(con))
  if (identical(orientation, "VECTORIZED")) {
    # convert to matrix and vector before writiing
    writeBin(as.vector(as.matrix(.data$signals)), con, size = 4)
  } else if (identical(orientation, "MULTIPLEXED")) {
    # transpose matrix before vectorizing
    writeBin(as.vector(t(as.matrix(.data$signals))), con, size = 4)
  }
}

write_vmrk <- function(.data,
                       filename) {
  vmrk_file <- paste0(filename,
                      ".vmrk")
  con <- file(vmrk_file,
              open = "w")
  on.exit(close(con))
  writeLines("Brain Vision Data Exchange Header File Version 2.0", con)
  writeLines(paste0("; Created using eegUtils http://craddm.github.io/eegUtils/"),
             con)
  writeLines("", con)

  new_markers <- list("Common Infos" = list(Codepage = "UTF-8",
                                            DataFile = paste0(filename, ".dat")),
                      "Marker Infos" = as.list(paste("Stimulus",
                                                     .data$events$event_type,
                                                     .data$events$event_onset,
                                                     1,
                                                     0,
                                                     sep = ",")))
  names(new_markers[["Marker Infos"]]) <-
    paste0("Mk", 1:nrow(.data$events))

  for (i in names(new_markers)) {
    writeLines(paste0("[", i,"]"), con)
    entries <-
      lapply(seq_along(new_markers[[i]]),
             function(x) paste0(names(new_markers[[i]][x]),
                                "=",
                                new_markers[[i]][x]))
    lapply(entries,
           writeLines,
           con)

    writeLines("", con)
  }
}
