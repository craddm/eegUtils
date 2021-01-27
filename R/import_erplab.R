#' Import from ERPLAB .erp files
#'
#' ERPLAB is a toolbox for Event Related Potential analysis associated with
#' EEGLAB. It exports files in Matlab v6.5 format, either with an .erp or .mat
#' extension. Note that ERPLAB files can contain data from multiple subjects,
#' individual subjects, or averaged over multiple subjects. The files do not
#' distinguish between these possibilities, so users will need to apply their
#' own knowledge of the contents to import them correctly.
#'
#' @param file_name Filename (and path if not in present working directory)
#' @param df_out Defaults to FALSE - outputs an object of class \code{eeg_epochs}. Set to
#'   TRUE to output a data frame.
#' @param participant_id By default, the filename will be used as the id of the
#'   participant.
#' @param recording By default, the filename will be used as the name of the
#'   recording.
#' @return An object of class \code{eeg_epochs} or a data.frame.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @examples
#' \dontrun{import_set("your_data.set")}
#' @export
import_erplab <- function(file_name,
                          df_out = FALSE,
                          participant_id = NULL,
                          recording = NULL) {

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Package \"R.matlab\" needed. Please install it.",
         call. = FALSE)
  }

  if (is.null(recording)) {
    recording <- basename(tools::file_path_sans_ext(file_name))
  }

  if (is.null(participant_id)) {
    participant_id <- basename(tools::file_path_sans_ext(file_name))
  }

  temp_dat <- R.matlab::readMat(file_name)
  times_data <- unlist(temp_dat$ERP["times", ,])
  n_times <- length(times_data)
  n_bins <- unlist(temp_dat$ERP["nbin", ,])
  n_chans <- unlist(temp_dat$ERP["nchan", ,])
  signals <- unlist(temp_dat$ERP["bindata", ,1])
  signals <- array(signals,
                   dim = c(n_chans,
                           n_times *
                           n_bins))
  signals <- t(signals)

  chan_info <- drop(Reduce(rbind, temp_dat$ERP["chanlocs",,]))
  chan_info <- apply(chan_info, 1, unlist)
  chan_info <- lapply(chan_info,
                      function(x) if (length(x) == 0) NA else x)
  chan_info <- tibble::as_tibble(chan_info)
  chan_info <- parse_chaninfo(chan_info)

  colnames(signals) <- chan_info$electrode
  srate <- unlist(temp_dat$ERP["srate", ,])

  timings <- tibble::tibble(time = rep(times_data, n_bins),
                            epoch = rep(1:n_bins, each = n_times),
                            sample = 1:nrow(signals))
  epochs <-
    tibble::new_tibble(list(epoch = 1:n_bins,
                            participant_id = rep(participant_id, n_bins),
                            recording = rep(recording, n_bins),
                            epoch_labels = unlist(temp_dat$ERP["bindescr", ,])),
                       nrow = n_bins,
                       class = "epoch_info")
  if (df_out) {
    return(dplyr::left_join(cbind(signals, timings), epochs))
  }

  eeg_epochs(tibble::as_tibble(signals),
             srate = srate,
             timings = timings,
             chan_info = chan_info,
             events = NULL,
             reference = NULL,
             epochs = epochs)
}
