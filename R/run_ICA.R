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
#' @param tol Convergence tolerance for fastica and infomax. Defaults to 1e-06.
#' @describeIn run_ICA Run ICA on an \code{eeg_epochs} object
#' @export

run_ICA.eeg_epochs <- function(data,
                               method = "sobi",
                               maxit = 500,
                               tol = 1e-6,
                               ...) {

  rank_check <- Matrix::rankMatrix(as.matrix(data$signals))
  #if (rank_check < n_channels) {
   # warning(paste("Data is rank deficient. Detected rank", rank_check))
    #pca <- prcomp(data$signals)
    #data$signals <- as.data.frame(pca$x[, 1:rank_check] %*% t(pca$rotation[, 1:rank_check]))
  #}

  if (method == "sobi") {
    if (!requireNamespace("JADE", quietly = TRUE)) {
      stop("Package \"JADE\" needed to use SOBI. Please install it.",
           call. = FALSE)
    }
    ica_obj <- sobi_ICA(data,
                        maxiter = maxit)
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
                              maxit = maxit,
                              tol = tol)
    } else if (method == "infomax") {
      ICA_out <- ica::icaimax(data$signals,
                              rank_check,
                              maxit = maxit,
                              fun = "ext",
                              tol = tol)
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
                       epochs = data$epochs,
                       algorithm = method)
    }
  ica_obj
}

#' SOBI ICA
#'
#' Internal function for running SOBI ICA on an \code{eeg_epochs} object
#'
#' @param data Data to be ICAed.
#' @param maxiter Maximum number of iterations of the joint diagonalization
#' @author A. Belouchrani and A. Cichocki. Adapted to R by Matt Craddock
#'   \email{matt@@mattcraddock.com}
#' @keywords internal

sobi_ICA <- function(data,
                     maxiter) {

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
  # # amp_matrix <- lapply(amp_matrix,
  # #                      function(x) scale(x, scale = FALSE))
   amp_matrix <- as.matrix(do.call("rbind", amp_matrix))
  #amp_matrix <- sweep(data$signals, 2, Matrix::colMeans(data$signals))
  SVD_amp <- svd(amp_matrix)

  ## get the psuedo-inverse of the diagonal matrix, multiply by singular
  ## vectors
   Q <- MASS::ginv(diag(SVD_amp$d), tol = 0) %*% t(SVD_amp$v) # whitening matrix
   amp_matrix <- Q %*% t(amp_matrix)
  #amp_matrix <- amp_matrix %*% SVD_amp$v

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
    M[, u:(u + n_channels - 1)] <- norm(Rxp, "F") * (Rxp)
  }

  epsil <- 1 / sqrt(N) / 100
  V <- JADE::rjd(t(M),
                 eps = epsil,
                 maxiter = maxiter)$V

  ## create mixing matrix for output
  mixing_matrix <- MASS::ginv(Q, tol = 0) %*% V

  unmixing_matrix <- MASS::ginv(mixing_matrix, tol = 0)
  # rescale vecs
  scaling <- sqrt(colMeans(mixing_matrix^2))

  unmixing_matrix <- sweep(unmixing_matrix, MARGIN = 1, scaling, `*`) # scaled weights

  mixing_matrix <- MASS::ginv(unmixing_matrix %*% diag(ncol(unmixing_matrix)),
                              tol = 0)

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

  ica_obj <- eeg_ICA(mixing_matrix = mixing_matrix,
                     unmixing_matrix = unmixing_matrix,
                     signals = S,
                     timings = data$timings,
                     events = data$events,
                     chan_info = data$chan_info,
                     srate = data$srate,
                     epochs = data$epochs,
                     algorithm = "sobi")
  ica_obj
}

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
    new_mixmat[, comps] <- 0
    source_sigs <- as.matrix(data$signals) %*%
      t(as.matrix(decomp$unmixing_matrix[, 1:ncomps]))

    new_dat <- as.matrix(new_mixmat) %*% t(source_sigs)
    new_dat <- as.data.frame(t(new_dat))
    names(new_dat) <- data$chan_info$electrode
    data$signals <- new_dat
    data
  }
}
