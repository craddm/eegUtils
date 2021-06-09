# Parse channel locations from HDF5 formatted EEGLAB
#'@keywords internal
parse_locshdf5 <- function(x) {
  chan_locs <- x[["chanlocs"]]
  fields <- names(chan_locs)
  file_struct <-
    lapply(fields,
           function(x) {
             chan_locs[[x]][1,]$dereference()
           })
  names(file_struct) <- fields

  data.frame(
    labels =
      vapply(file_struct$labels,
             function(x)
               intToUtf8(x[1, ]),
             character(1)),
    radius =
      vapply(file_struct$radius,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    ref =
      vapply(file_struct$ref,
             function(x)
               intToUtf8(x[1, ]),
             character(1)),
    sph.phi =
      vapply(file_struct$sph_phi,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    sph.radius =
      vapply(file_struct$sph_radius,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    sph.theta =
      vapply(file_struct$sph_theta,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    theta =
      vapply(file_struct$theta,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    type =
      vapply(file_struct$type,
             function(x)
               intToUtf8(x[]),
             character(1)),
    urchan =
      vapply(file_struct$urchan,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    X =
      vapply(file_struct$X,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    Y =
      vapply(file_struct$Y,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1)),
    Z =
      vapply(file_struct$Z,
             function(x) {
               if (length(x$dims) == 1)
                 NA
               else
                 x[1, ]
             },
             numeric(1))
  )
}

parse_events_hdf5 <- function(x) {
  get_events <- x[["event"]]
  fields <- names(get_events)
  file_struct <-
    lapply(fields,
           function(x) {
             get_events[[x]][1,]$dereference()
           })
  names(file_struct) <- fields
  event_table <-
    as.data.frame(lapply(file_struct,
                         function(x)
                           unlist(lapply(x,
                                         function(z)
                                           z$read()))))
  tibble::as_tibble(event_table)
}

read_hdf5_set <- function(x,
                          recording,
                          participant_id,
                          drop_custom = TRUE) {
  eeg_set <- hdf5r::H5File$new(x,
                               mode = "r")
  signals <- eeg_set[["data"]]$read()
  dim_signals <- dim(signals)
  dim(signals) <- c(dim_signals[1], dim_signals[2] * dim_signals[3])

  chaninfo <- parse_chaninfo(parse_locshdf5(eeg_set))
  signals <- t(signals)
  colnames(signals) <- chaninfo$electrode
  signals <- tibble::as_tibble(signals)

  srate <- eeg_set[["srate"]]$read()
  reference <- intToUtf8(eeg_set[["ref"]]$read())
  times <- eeg_set[["times"]]$read()
  n_epochs <- eeg_set[["trials"]]$read()

  timings <- tibble::tibble(
    time = rep(times / 1000, n_epochs),
    epoch = rep(1:n_epochs, each = length(times)),
    sample = 1:nrow(signals)
  )

  event_table <- parse_events_hdf5(eeg_set)

  event_table$event_time <- (event_table$latency - 1) / srate

  if (any(event_table$latency %% 1 > 0)) {
    message("Rounding non-integer event sample latencies...")
    event_table$latency <- round(event_table$latency)
  }

  # EEGLAB stores latencies in samples starting from 1, my event_time is in
  # seconds, starting from 0
  event_table$event_time <- (event_table$latency - 1) / srate

  std_cols <- c("latency",
                "event_time",
                "type",
                "epoch")

  if (drop_custom & any(!colnames(event_table) %in% std_cols)) {
    message("Dropping custom columns...")
    event_table <- event_table[, std_cols]
  }

  col_check <-  colnames(event_table) %in% c("event_type", "event_onset")

  if (any(col_check)) {
    dupe_checks <- colnames(event_table)[col_check]
    dupe_labs <- paste0(dupe_checks, ".x")
  }

  #need to build in check for duplicate columns
  event_table <- dplyr::rename(event_table,
                               event_type = "type",
                               event_onset = "latency")

  event_table$time <- NA
  event_table$time <- timings[which(timings$sample %in% event_table$event_onset,
                                    arr.ind = TRUE), ]$time

  eeg_set$close_all()

  epochs <-
    tibble::new_tibble(
      list(
        epoch = 1:n_epochs,
        participant_id = rep(participant_id, n_epochs),
        recording = rep(recording, n_epochs)
      ),
      nrow = n_epochs,
      class = "epoch_info"
    )

  eeg_epochs(
    data = signals,
    srate = srate,
    chan_info = chaninfo,
    reference = list(ref_chans = reference,
                     excluded = NULL),
    timings = timings,
    epochs = epochs,
    events = event_table
  )
}
