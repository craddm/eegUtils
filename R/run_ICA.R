#' ICA for EEG data
#'
#' Implements SOBI ICA. Although SOBI exists in R already, thanks to the JADE package, it doesn't respect epoch boundaries when computing correlations across different time lags. This is a port of SOBI from EEGLAB. Currently only works on epoched data.
#'
#' @param df Data frame to be ICAed.
#' @param method Only SOBI is currently implemented, so this is ignored.
#'
#' @import dplyr
#' @import JADE
#' @export

run_ICA <- function (df, method = "sobi") {

  n_epochs <- length(unique(df$epoch))
  n_channels <- length(unique(df$electrode))
  n_times <- length(unique(df$time))

  ##number of lags at which to assess autocovariance matrices
  n_lags <- min(100,ceiling(n_times/3))

  ## reshape amplitude X electrode to a square matrix before doing SVD
  df <- df %>% spread(electrode, amplitude)

  ## Pre-whiten the data using the SVD.
  ## zero-mean columns and get SVD. NB: should probably edit this to zero mean *epochs*
  amp_matrix <- scale(data.matrix(df[ , (ncol(df) - n_channels + 1):ncol(df)]), scale = FALSE)
  SVD_amp <- svd(amp_matrix)

  ## get the psuedo-inverse of the diagonal matrix, multiply by right singular vectors
  Q <- MASS::ginv(diag(SVD_amp$d)) %*% t(SVD_amp$v)
  amp_matrix <- Q %*% t(amp_matrix)

  ## reshape to reflect epoching structure
  dim(amp_matrix) <- c(n_channels, n_times, n_epochs)

  ## vectorise this loop ASAP.
  k <- 1
  pm <- n_channels * n_lags
  N <- n_times
  M <- matrix(NA,nrow = n_channels, ncol = pm)
  for (u in seq(1, pm, n_channels)) {
    k <- k + 1
    for (tlag in 1:n_epochs) {
      if (tlag == 1) {
        Rxp = amp_matrix[ , k:N, tlag] %*% t(amp_matrix[ , 1:(N - k + 1), tlag]) / (N - k + 1) / n_epochs }
      else {
        Rxp = Rxp + amp_matrix[, k:N, tlag] %*% t(amp_matrix[ , 1:(N - k + 1), tlag]) / (N - k + 1) / n_epochs
      }
      M[ , u:(u + n_channels - 1)] <- norm(Rxp, 'F') * Rxp #  % Frobenius norm =
    }
  }

  dim(M) <- c(n_channels, n_channels, pm/n_channels)

  ## do joint diagonalization of matrices using rjd from JADE
  M_rjd <- JADE::rjd(M, maxiter = 300, eps = 1 / sqrt(N) / 100)

  ## create mixing matrix for output
  mixing_matrix <- data.frame(MASS::ginv(Q) %*% M_rjd$V)
  names(mixing_matrix) <- 1:64
  mixing_matrix$electrode <- names(df[ , (ncol(df) - n_channels + 1):ncol(df)])
  dim(amp_matrix) <- c(n_channels, n_times*n_epochs)
  S <- t(M_rjd$V) %*% amp_matrix
  return(list("mixing_matrix" = mixing_matrix, "comp_activations" = t(S)))
}
