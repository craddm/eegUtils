#' Generalized eigenvalue decomposition based methods for EEG data
#'
#' Implements a selection of Generalized Eigenvalue based decomposition methods
#' for EEG signals. Intended for isolating oscillations at specified frequencis,
#' decomposing channel-based data into distinct components reflecting distinct
#' or combinations of sources of oscillatory signals. Currently only the
#' spatio-spectral decomposition method (Nikulin et al, 2011) is implemented.
#'
#' @param data An \code{eeg_data} object
#' @param ... Additional parameters
#' @references Cohen, M. X. (2016). Comparison of linear spatial filters for
#'   identifying oscillatory activity in multichannel data. BioRxiv, 097402.
#'   https://doi.org/10.1101/097402
#'
#'   Haufe, S., Dähne, S., & Nikulin, V. V.
#'   (2014). Dimensionality reduction for the analysis of brain oscillations.
#'   NeuroImage, 101, 583–597. https://doi.org/10.1016/j.neuroimage.2014.06.073
#'
#'   Nikulin, V. V., Nolte, G., & Curio, G. (2011). A novel method for reliable
#'   and fast extraction of neuronal EEG/MEG oscillations on the basis of
#'   spatio-spectral decomposition. NeuroImage, 55(4), 1528–1535.
#'   https://doi.org/10.1016/j.neuroimage.2011.01.057

eeg_decomp <- function(data, ...) {
  UseMethod("eeg_decomp", data)
}

eeg_decomp.default <- function(data, ...) {
  stop("Not implemented for objects of class ", class(data))
}

#' @param sig_range Vector with two inputs, the lower and upper bounds of the frequency range of interest
#' @param noise_range Range of frequencies to be considered noise (e.g. bounds of flanker frequencies)
#' @param method Type of decomposition to apply. Currently only "ssd" is supported.
#' @describeIn eeg_decomp method for \code{eeg_epochs} objects
eeg_decomp.eeg_epochs <- function(data,
                                  sig_range,
                                  noise_range = NULL,
                                  method = "ssd") {


  data <- switch(method,
                 "ssd" = run_SSD(data,
                                 sig_range,
                                 noise_range),
                 "ress" = run_SSD(data,
                                  sig_range,
                                  noise_range,
                                  RESS = TRUE))
  class(data) <- c("eeg_ICA", "eeg_epochs")
  data
}

#' Internal function for running SSD algorithm
#'
#' @param data \code{eeg_epochs} object to be decomposed
#' @param sig_range Frequency range of the signal of interest
#' @param noise_range Frequency range of the noise
#' @keywords internal

run_SSD <- function(data,
                    sig_range,
                    noise_range,
                    RESS = FALSE) {

  if (!requireNamespace("geigen", quietly = TRUE)) {
    stop("Package \"geigen\" needed for SSD. Please install it.",
         call. = FALSE)
  }

  signal <- iir_filt(data,
                     low_freq = sig_range[1],
                     high_freq = sig_range[2],
                     filter_order = 2)
  noise <- iir_filt(data,
                    low_freq = noise_range[1],
                    high_freq = noise_range[2],
                    filter_order = 2)
  noise <- iir_filt(data,
                    low_freq = (sig_range[2] + noise_range[2]) / 2,
                    high_freq = (sig_range[1] + noise_range[1]) / 2,
                    filter_order = 2)
  # Calculate covariance respecting the epoching structure of the data
  cov_sig <- cov_epochs(signal)
  cov_noise <- cov_epochs(noise)

  eig_sigs <- base::eigen(cov_sig)
  # Get the rank of the covariance matrix and select only as many components as
  # there are ranks
  rank_sig <- Matrix::rankMatrix(cov_sig)

  if (rank_sig < ncol(cov_sig)) {
    message("Input data is not full rank; returning ",
            rank_sig,
            "components")
  }

  M <- eig_sigs$vectors[, 1:rank_sig] %*% (diag(eig_sigs$values[1:rank_sig] ^ -0.5))

  C_s_r <- t(M) %*% cov_sig %*% M
  C_n_r <- t(M) %*% cov_noise %*% M
  ged_v <- geigen::geigen(C_s_r, C_s_r + C_n_r) # this one needs to be sorted
  lambda <- sort(ged_v$values, decreasing = TRUE)
  W <- ged_v$vectors[,order(ged_v$values, decreasing = TRUE)]
  W <- M %*% W
  data$mixing_matrix <- (cov_sig %*% W) %*% solve(t(W) %*% cov_sig %*% W)
  data$mixing_matrix <- as.data.frame(data$mixing_matrix)
  names(data$mixing_matrix) <- paste0("Comp", 1:ncol(data$mixing_matrix))
  data$mixing_matrix$electrode <- names(data$signals)

  data$unmixing_matrix <- as.data.frame(W)
  names(data$unmixing_matrix) <- paste0("Comp", 1:ncol(W))
  data$unmixing_matrix$electrode <- data$mixing_matrix$electrode

  if (RESS) {
    data$signals <- as.data.frame(as.matrix(data$signals) %*% W)
    names(data$signals) <- paste0("Comp", 1:ncol(W))
    return(data)
  }
  data$signals <- as.data.frame(as.matrix(signal$signals) %*% W)
  names(data$signals) <- paste0("Comp", 1:ncol(W))
  data

}

#' Covariance of epoched data
#'
#' Calculate covariance of each epoch, then average
#'
#' @param data epoched data to calculate covariance for
#' @keywords internal

cov_epochs <- function(data) {

  if (!is.eeg_epochs(data)) {
    stop("This is not an eeg_objects object.")
  }
  # Data is converted to a 3D matrix (n_times X n_epochs X n_channels),
  # covariance is calculated for each epoch then averaged over
  full_cov <- rowMeans(apply(conv_to_mat(data),
                             2,
                             stats::cov))
  dim(full_cov) <- c(ncol(data$signals),
                     ncol(data$signals))
  as.matrix(full_cov)
}
