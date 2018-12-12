#' Independent Component Analysis for EEG data
#'
#' Performs Independent Component Analysis for EEG data. Currently only
#' available with on epoched data. Implements three different methods of ICA -
#' fastica, extended Infomax, and Second-Order Blind Identification (SOBI).
#'
#' @param data Data to be ICAed.
#' @param ... Other parameters passed to function.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @return An \code{eeg_ICA} object containing an ICA decomposition
#' @importFrom MASS ginv
#' @export

run_ICA <- function(data, ...) {
  UseMethod("run_ICA", data)
}

#' @param method "sobi" (default), "fastica", or "infomax".
#' @param maxit Maximum number of iterations of the Infomax and Fastica ICA
#'   algorithms.
#' @describeIn run_ICA Run ICA on an \code{eeg_epochs} object
#' @export

run_ICA.eeg_epochs <- function(data,
                               method = "sobi",
                               maxit = 500,
                               ...) {

  rank_check <- Matrix::rankMatrix(as.matrix(data$signals))
  #if (rank_check < n_channels) {
   # warning(paste("Data is rank deficient. Detected rank", rank_check))
    #pca <- prcomp(data$signals)
    #data$signals <- as.data.frame(pca$x[, 1:rank_check] %*% t(pca$rotation[, 1:rank_check]))
  #}

  if (method == "sobi") {
    ica_obj <- sobi_ICA(data)
  } else {

    if (!requireNamespace("ica", quietly = TRUE)) {
      stop("Package \"ica\" needed to use infomax or fastica. Please install it.",
           call. = FALSE)
    }

    if (rank_check < ncol(data$signals)) {
      warning(paste("Data is rank deficient. Detected rank", rank_check))
    }

    if (method == "fastica") {
      ICA_out <- ica::icafast(data$signals,
                              rank_check,
                              maxit = maxit)
    } else if (method == "infomax") {
      ICA_out <- ica::icaimax(data$signals,
                              rank_check,
                              maxit = maxit,
                              fun = "ext")
    } else {
      stop("Unknown method; available methods are sobi, fastica, and infomax.")
    }

    ICA_out$S <- as.data.frame(ICA_out$S)
    names(ICA_out$S) <- paste0("Comp", 1:ncol(ICA_out$S))
    mixing_matrix <- as.data.frame(ICA_out$M)
    names(mixing_matrix) <- paste0("Comp", 1:ncol(ICA_out$S))
    mixing_matrix$electrode <- names(data$signals)
    unmixing_matrix <- as.data.frame(t(ICA_out$W))
    names(unmixing_matrix) <- paste0("Comp", 1:ncol(ICA_out$S))
    unmixing_matrix$electrode <- names(data$signals)

    ica_obj <- eeg_ICA(mixing_matrix = mixing_matrix,
                    unmixing_matrix = unmixing_matrix,
                    signals = ICA_out$S,
                    timings = data$timings,
                    events = data$events,
                    chan_info = data$chan_info,
                    srate = data$srate,
                    continuous = FALSE)
  }
  ica_obj
}

#' SOBI ICA
#'
#' Internal function for running SOBI ICA on an \code{eeg_epochs} object
#'
#' @param data Data to be ICAed.
#' @keywords internal

sobi_ICA <- function(data) {

  n_epochs <- length(unique(data$timings$epoch))
  n_channels <- ncol(data$signals)
  n_times <- length(unique(data$timings$time))

  ##number of lags at which to assess autocovariance matrices
  ## 100 is the default; only switches to smaller n if epochs are too short
  n_lags <- min(100, ceiling(n_times / 3))

  ## Pre-whiten the data using the SVD. zero-mean columns and get SVD. NB:
  ## should probably edit this to zero mean *epochs*
  #do by epochs
  amp_matrix <- split(data$signals, data$timings$epoch)
  amp_matrix <- lapply(amp_matrix,
                       function(x) sweep(x,
                                         2,
                                         Matrix::colMeans(x))
                       )
  amp_matrix <- as.matrix(do.call("rbind", amp_matrix))
  SVD_amp <- svd(amp_matrix)

  ## get the psuedo-inverse of the diagonal matrix, multiply by singular
  ## vectors
  Q <- MASS::ginv(diag(SVD_amp$d), tol = 0) %*% t(SVD_amp$v) # whitening matrix
  amp_matrix <- Q %*% t(amp_matrix)

  ## reshape to reflect epoching structure
  dim(amp_matrix) <- c(n_channels,
                       n_times,
                       n_epochs)

  ## vectorise this loop if possible?
  k <- 1
  pm <- n_channels * n_lags
  N <- n_times
  M <- matrix(NA, nrow = n_channels, ncol = pm)

  tmp_fun <- function(amps,
                      k,
                      N) {
    amps[, k:N] %*% t(amps[, 1:(N - k + 1)]) / (N - k + 1)
  }

  for (u in seq(1, pm, n_channels)) {
    k <- k + 1
    Rxp <- lapply(seq(1, n_epochs),
                  function(x) tmp_fun(amp_matrix[, , x],
                                      k,
                                      N))
    Rxp <- Reduce("+",
                  Rxp)
    Rxp <- Rxp / n_epochs
    M[, u:(u + n_channels - 1)] <- norm(Rxp, "F") * (Rxp) #  Frobenius norm
  }

  epsil <- 1 / sqrt(N) / 100
  V <- joint_diag(M,
                  epsil,
                  n_channels,
                  pm) # V = sphered unmixing matrix

  ## create mixing matrix for output
  mixing_matrix <- MASS::ginv(Q, tol = 0) %*% V
  unmixing_matrix <- t(V)

  dim(amp_matrix) <- c(n_channels,
                       n_times * n_epochs)

  S <- unmixing_matrix %*% t(as.matrix(data$signals))
  S <- as.data.frame(t(S))
  names(S) <- paste0("Comp", 1:ncol(S))

  mixing_matrix <- as.data.frame(mixing_matrix)
  names(mixing_matrix) <- paste0("Comp", 1:ncol(S))
  mixing_matrix$electrode <- names(data$signals)
  unmixing_matrix <- as.data.frame(unmixing_matrix)
  names(unmixing_matrix) <- paste0("Comp", 1:ncol(S))
  unmixing_matrix$electrode <- names(data$signals)

  ica_obj <- list("mixing_matrix" = mixing_matrix,
                  "unmixing_matrix" = unmixing_matrix,
                  "signals" = S,
                  "timings" = data$timings,
                  "events" = data$events,
                  "chan_info" = data$chan_info,
                  "srate" = data$srate,
                  "continuous" = FALSE)
  class(ica_obj) <- c("eeg_ICA", "eeg_epochs")
  return(ica_obj)
}

#' Joint diagonalization for SOBI
#'
#' Drop-in replacement for JADE::rjd, which often fails to converge
#'
#' @param M Numeric matrix
#' @param eps Convergence tolerance
#' @param pm number of channels * number of timepoints
#' @param n_channels number of channels
#' @keywords internal

joint_diag <- function(M,
                       eps,
                       n_channels,
                       pm) {

  epsil <- eps
  continue <- TRUE
  V <- diag(n_channels)
  step_n <- 0

  while (isTRUE(continue)) {
    continue <- FALSE
    for (p in 1:(n_channels - 1)) {
      for (q in ((p + 1):n_channels)) {

        P_seq <- seq(p, pm, n_channels)
        Q_seq <- seq(q, pm, n_channels)

        # Perform Givens rotation
        g <- rbind(M[p, P_seq] - M[q, Q_seq],
                   M[p, Q_seq] + M[q, P_seq],
                   1i * (M[q, P_seq] - M[p, Q_seq]))
        eigs <- eigen(Re(tcrossprod(g))) # %*% t(g)))
        sort_eigs <- sort(eigs$`values`, index.return = T)
        angles <- eigs$vectors[, sort_eigs$ix[3]]
        angles <- sign(angles[1]) * angles
        c_r <- sqrt(0.5 + angles[1] / 2)
        sr <- Re(0.5 * (angles[2] - 1i * angles[3]) / c_r)
        sc <- Conj(sr)
        conv_check <- abs(sr) > epsil
        continue <- (continue | conv_check)
        if (conv_check) {# Update the M and V matrices
          colp <- M[, P_seq]
          colq <- M[, Q_seq]
          M[, P_seq] <- c_r * colp + sr * colq
          M[, Q_seq] <- c_r * colq - sc * colp
          rowp <- M[p, ]
          rowq <- M[q, ]
          M[p, ] <- c_r * rowp + sc * rowq
          M[q, ] <- c_r * rowq - sr * rowp
          temp <- V[, p]
          V[, p] <- c_r * V[, p] + sr * V[, q]
          V[, q] <- c_r * V[, q] - sc * temp
        }
      }
    }
    step_n <- step_n + 1
    message("Iteration = ", step_n, "\n")
  }
  V
}


#

apply_ica <- function(data, ...) {
  UseMethod("apply_ica", data)
}

apply_ica.eeg_ICA <- function(data, comps, ...) {
  ncomps <- ncol(data$mixing_matrix)
  new_mixmat <- data$mixing_matrix[1:(ncomps - 1)]
  new_mixmat[ , comps] <- 0
  new_dat <- as.matrix(new_mixmat) %*% t(as.matrix(data$signals))
  new_dat <- as.data.frame(t(new_dat))
  names(new_dat) <- data$chan_info$electrode
  out <-
    eeg_data(new_dat,
           data$srate,
           data$events,
           data$chan_info,
           data$timings,
           data$continuous,
           reference = NULL
          )
  class(out) <- c("eeg_epochs", "eeg_data")
  out
}

apply_ica.eeg_epochs <- function(data, decomp, comps, ...) {
  if (!is.eeg_ICA(decomp)) {
    stop("An eeg_ICA object must be supplied as decomp.")
  } else {
    ncomps <- ncol(decomp$mixing_matrix) - 1
    new_mixmat <- decomp$mixing_matrix[1:(ncomps)]
    new_mixmat[ , comps] <- 0
    source_sigs <- as.matrix(data$signals) %*%
      as.matrix(decomp$unmixing_matrix[, 1:ncomps])

    new_dat <- as.matrix(new_mixmat) %*% t(source_sigs)
    new_dat <- as.data.frame(t(new_dat))
    names(new_dat) <- data$chan_info$electrode
    data$signals <- new_dat
    data
  }
}
