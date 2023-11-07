ft_raw <- function(data,
                   struct_names,
                   participant_id,
                   recording,
                   verbose) {

  if (verbose) message("Importing ft_raw (epoched) dataset.")

  sig_names <- unlist(data["label", , ],
                      use.names = FALSE)

  times <- unlist(data["time", , ],
                  use.names = FALSE)

  if ("fsample" %in% struct_names) {
    srate <- unlist(data["fsample", , ],
                    use.names = FALSE)
  } else {
    srate <- NULL
  }

  n_trials <- length(data["trial", , ][[1]])
  signals <- purrr::flatten(purrr::flatten(data["trial", , ]))
  signals <- tibble::as_tibble(t(do.call(cbind, signals)))
  names(signals) <- sig_names
  timings <- tibble::tibble(epoch = rep(seq(1, n_trials),
                                        each = length(unique(times))),
                            time = times)

  epochs <- tibble::tibble(epoch = 1:n_trials,
                           participant_id = participant_id,
                           recording = recording)

  if ("elec" %in% struct_names) {
    if (verbose) message("Importing channel information.")
    chan_info <- ft_chan_info(
      unlist(data["elec", , ],
             recursive = FALSE)
    )
    chan_info$electrode <- sig_names
    chan_info <- validate_channels(chan_info)
  } else if ("grad" %in% struct_names) {
    message("MEG channel locations found; import of locations not currently supported.")

    chan_info <- data.frame(electrode = sig_names)
    chan_info <- validate_channels(chan_info)
  } else {
    if (verbose) message("No EEG channel information found.")
    chan_info <- NULL
  }

  eeg_epochs(data = signals,
             srate = srate,
             timings = timings,
             chan_info = chan_info,
             epochs = epochs)
}



proc_dimord <- function(dimord) {
  dimord <- unlist(strsplit(dimord, "_"))
  dimord <- gsub("rpt", "epoch", dimord)
  dimord <- gsub("chan", "electrode", dimord)
  dimord <- gsub("freq", "frequency", dimord)
  dimord
}
ft_comp <- function(data) {
  mixing_matrix <- unlist(data["topo", , ],
                          use.names = FALSE)
  unmixing_matrix <- unlist(data["unmixing", , ],
                            use.names = FALSE)
  chan_names <- unlist(data["topolabel", , ],
                       use.names = FALSE)
}

ft_freq <- function(data,
                    dimensions,
                    struct_names,
                    participant_id,
                    recording,
                    verbose) {

  if (verbose) message("Importing ft_freq dataset.")

  frequencies <- as.vector(unlist(data["freq", , ]))

  if (!("powspctrm" %in% struct_names)) {
    stop("Only import of power data currently supported from ft_freq data.")
  }

  if ("fsample" %in% struct_names) {
    srate <- unlist(data["fsample", ,],
                    use.names = FALSE)
  } else {
    srate <- NULL
  }

  signals <- data[[which(struct_names == "powspctrm")]]

  sig_names <- unlist(data["label", , ],
                      use.names = FALSE)

  if ("time" %in% dimensions) {
    times <- unlist(data["time", ,],
                    use.names = FALSE)
  }

  freq_info <- list(freqs = frequencies,
                    method = "unknown",
                    output = "power",
                    baseline = "none")

  if (all(c("epoch", "time") %in% dimensions)) {
    n_trials <- dim(signals)[[1]]
    dimnames(signals) <- list(epoch = 1:n_trials,
                              electrode = sig_names,
                              frequency = frequencies,
                              time = times)
  } else {
    n_trials <- 1
    dimnames(signals) <- list(electrode = sig_names,
                              frequency = frequencies,
                              time = times)
  }

  epochs <- tibble::tibble(epoch = 1:n_trials,
                           participant_id = participant_id,
                           recording = recording)

  timings <- tibble::tibble(epoch = rep(seq(1, n_trials),
                                        each = length(unique(times))),
                            time = rep(times,
                                       n_trials))

  if ("elec" %in% struct_names) {
    if (verbose) message("Importing channel information.")
    chan_info <- ft_chan_info(
      unlist(data["elec", ,],
             recursive = FALSE)
    )
    chan_info$electrode <- sig_names
    chan_info <- validate_channels(chan_info)
  } else {
    if (verbose) message("No EEG channel information found.")
    chan_info <- NULL
  }

  eeg_tfr(data = signals,
          dimensions = dimensions,
          srate = srate,
          events = NULL,
          chan_info = chan_info,
          reference = NULL,
          timings = timings,
          epochs = epochs,
          freq_info = freq_info)
}


ft_chan_info <- function(data) {
  # chan_info <- unlist(data["elec", ,],
  #                     recursive = FALSE)
  chan_info <- data[[1]]
  colnames(chan_info) <- c("cart_x",
                           "cart_y",
                           "cart_z")
  chan_info <- tibble::as_tibble(chan_info)
  sph_coords <- cart_to_spherical(chan_info)
  chan_info <- cbind(chan_info,
                     sph_coords)
  chan_info[, c("x", "y")] <- project_elecs(chan_info)
  chan_info
}
