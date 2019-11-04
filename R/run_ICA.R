#' Independent Component Analysis for EEG data
#'
#' Performs Independent Component Analysis for EEG data. Currently only
#' available with on epoched data. Implements three different methods of ICA -
#' fastica, extended Infomax, and Second-Order Blind Identification (SOBI).
#'
#' @section Notes on ICA usage:
#'
#'   It is recommended to mean-centre your data appropriately before running
#'   ICA. The implementations of FASTICA and extended-Infomax from the `ica`
#'   package, and of SOBI ICA have this as an option which is enabled by
#'   default, while the implementation of FASTICA in the fICA package enforces
#'   mean-centring of the columns of the data. With epoched data, it is
#'   recommended to centre each epoch on zero, rather than centre on the overall
#'   channel mean. This can be achieved with the `rm_baseline()` function. SOBI
#'   ICA will do this automatically, whereas the other ICA implementations will
#'   centre on the channel means, not the epoch means.
#'
#' @param data Data to be ICAed.
#' @param ... Other parameters passed to function.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @return An \code{eeg_ICA} object containing an ICA decomposition
#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix
#' @examples
#' run_ICA(demo_epochs)
#' run_ICA(demo_epochs, pca = 10)
#' @export

run_ICA <- function(data, ...) {
  UseMethod("run_ICA", data)
}

#' @param method "sobi" (default), "fastica", or "infomax".
#' @param maxit Maximum number of iterations of the Infomax and Fastica ICA
#'   algorithms.
#' @param tol Convergence tolerance for fastica and infomax. Defaults to 1e-06.
#' @param pca Reduce the number of dimensions using PCA before running ICA.
#'   Numeric,  >1 and < number of channels
#' @param centre Defaults to TRUE. Centre the data on zero by subtracting the
#'   column mean. See notes on usage.
#' @param alg Use "gradient" descent or "newton" algorithm for extended infomax.
#'   Defaults to "gradient". Ignored if method != "infomax".
#' @param rateanneal Annealing rate for extended infomax. Ignored if method !=
#'   "infomax".
#' @param rate Learning rate for extended infomax. Ignored if method !=
#'   "infomax".
#' @describeIn run_ICA Run ICA on an \code{eeg_epochs} object
#' @importFrom stats cov
#' @export

run_ICA.eeg_epochs <- function(data,
                               method = "sobi",
                               maxit = 1000,
                               tol = 1e-6,
                               pca = NULL,
                               centre = TRUE,
                               alg = "gradient",
                               rateanneal = c(60, .9),
                               rate = 0.1,
                               ...) {

  if (!is.null(pca)) {
    message("Reducing data to ", pca,
            " dimensions using PCA.")
    orig_chans <- channel_names(data)
    pca_decomp <- eigen(stats::cov(data$signals))$vectors
    data$signals <- as.data.frame(as.matrix(data$signals) %*% pca_decomp[, 1:pca])
    pca_flag <- TRUE
  } else {
    pca <- ncol(data$signals)
    pca_flag <- FALSE
  }

  rank_check <- Matrix::rankMatrix(as.matrix(data$signals))

  if (rank_check < ncol(data$signals)) {
    warning(paste("Data is rank deficient. Detected rank", rank_check))
  }

  if (!(method %in% c("sobi", "fastica", "infomax", "fica"))) {
    stop("Unknown method; available methods are sobi, fastica, infomax, and fica.")
  }

  if (method == "fica") {
    if (!requireNamespace("fICA", quietly = TRUE)) {
      stop("Package \"fICA\" needed to use the fICA implementation of fastica. Please install it.",
           call. = FALSE)
    }
    message("Running fastica (fICA).")

    ICA_out <- fICA::fICA(as.matrix(data$signals),
                          eps = tol,
                          maxiter = maxit,
                          method = "sym2")

    colnames(ICA_out$S) <- sprintf("Comp%03d", 1:pca)
    ICA_out$S <- tibble::as_tibble(ICA_out$S[, 1:pca])

    ICA_out$W <- ICA_out$W[, 1:pca]
    mixing_matrix <- MASS::ginv(ICA_out$W, tol = 0)
    mixing_matrix <- pca_decomp[, 1:pca] %*% mixing_matrix

    unmixing_matrix <- as.data.frame(MASS::ginv(mixing_matrix, tol = 0))
    mixing_matrix <- as.data.frame(mixing_matrix)

    names(mixing_matrix) <- sprintf("Comp%03d", 1:pca)
    mixing_matrix$electrode <- orig_chans
    names(unmixing_matrix) <- orig_chans
    unmixing_matrix$Component <- sprintf("Comp%03d", 1:pca)
  } else {

    if (method == "sobi") {
      if (!requireNamespace("JADE", quietly = TRUE)) {
        stop("Package \"JADE\" needed to use SOBI. Please install it.",
             call. = FALSE)
      }

      message("Running SOBI ICA.")
      ICA_out <- sobi_ICA(data,
                          maxiter = maxit,
                          tol = tol,
                          pca = pca,
                          centre = centre)
    } else {

      if (!requireNamespace("ica", quietly = TRUE)) {
        stop("Package \"ica\" needed to use infomax or fastica. Please install it.",
             call. = FALSE)
      }

      if (method == "fastica") {
        message("Running fastica (ica).")
        ICA_out <- ica::icafast(data$signals,
                                nc = rank_check,
                                maxit = maxit,
                                tol = tol,
                                center = centre)
      } else if (method == "infomax") {
        message("Running extended-Infomax (ica).")
        ICA_out <- ica::icaimax(data$signals,
                                nc = rank_check,
                                maxit = maxit,
                                fun = "ext",
                                tol = tol,
                                center = centre,
                                alg = alg,
                                rateanneal = rateanneal,
                                rate = rate)
      }
    }

      ICA_out$S <- as.data.frame(ICA_out$S)
      names(ICA_out$S) <- sprintf("Comp%03d", 1:ncol(ICA_out$S))

      if (pca_flag) {
        mixing_matrix <- as.data.frame(pca_decomp[, 1:pca] %*% ICA_out$M)
        unmixing_matrix <- as.data.frame(MASS::ginv(as.matrix(mixing_matrix), tol = 0))
        names(mixing_matrix) <- sprintf("Comp%03d", 1:pca)
        names(unmixing_matrix) <- orig_chans
        mixing_matrix$electrode <- orig_chans
        unmixing_matrix$Component <- sprintf("Comp%03d", 1:pca)
      } else {
        mixing_matrix <- as.data.frame(ICA_out$M)
        names(mixing_matrix) <- sprintf("Comp%03d", 1:ncol(ICA_out$M))
        mixing_matrix$electrode <- names(data$signals)
        unmixing_matrix <- as.data.frame(ICA_out$W)
        names(unmixing_matrix) <- names(data$signals)
        unmixing_matrix$Component <- sprintf("Comp%03d", 1:ncol(ICA_out$S))
      }
  }

  ica_obj <- eeg_ICA(mixing_matrix = mixing_matrix,
                     unmixing_matrix = unmixing_matrix,
                     signals = ICA_out$S,
                     timings = data$timings,
                     events = data$events,
                     chan_info = data$chan_info,
                     srate = data$srate,
                     epochs = data$epochs,
                     algorithm = method)
  ica_obj
}

#' SOBI ICA
#'
#' Internal function for running SOBI ICA on an \code{eeg_epochs} object
#'
#' @param data Data to be ICAed.
#' @param maxiter Maximum number of iterations of the joint diagonalization
#' @param tol convergence tolerance.
#' @param pca Number of PCA components.
#' @param centre Mean centre signals.
#' @param verbose Print informative messages.
#' @author A. Belouchrani and A. Cichocki. Adapted to R by Matt Craddock
#'   \email{matt@@mattcraddock.com}
#' @keywords internal

sobi_ICA <- function(data,
                     maxiter,
                     tol,
                     pca,
                     centre,
                     verbose = TRUE) {

  #n_epochs <- length(unique(data$timings$epoch))
  n_epochs <- nrow(epochs(data))
  n_channels <- ncol(data$signals)
  n_times <- length(unique(data$timings$time))

  ##number of lags at which to assess autocovariance matrices
  ## 100 is the default; only switches to smaller n if epochs are too short
  n_lags <- min(100,
                ceiling(n_times / 3))

  ## Pre-whiten the data using the SVD. zero-mean columns and get SVD. NB:
  ## should probably edit this to zero mean *epochs*
  #do by epochs

  if (centre){
    # centre the data on zero.
    data <- rm_baseline(data,
                        verbose = FALSE)
  }

  SVD_amp <- svd(data$signals)

  ## get the psuedo-inverse of the diagonal matrix, multiply by singular
  ## vectors
  Q <- tcrossprod(MASS::ginv(diag(SVD_amp$d),
                             tol = 0),
                  SVD_amp$v) # whitening matrix
  amp_matrix <- tcrossprod(Q,
                           as.matrix(data$signals))

  ## reshape to reflect epoching structure
  dim(amp_matrix) <- c(n_channels,
                       n_times,
                       n_epochs)

  ## vectorise this loop if possible?
  k <- -1
  pm <- n_channels * n_lags
  N <- n_times
  M <- matrix(NA,
              nrow = n_channels,
              ncol = pm)

  for (u in seq(1, pm, n_channels)) {
    k <- k + 1
    # Rxp <- lapply(seq(1, n_epochs),
    #               function(x) tmp_fun(amp_matrix[, , x],
    #                                   k,
    #                                   N))
    # Rxp <- Reduce("+",
    #               Rxp)
    # Rxp <- Rxp / n_epochs
    # M[, u:(u + n_channels - 1)] <- norm(Rxp, "F") * (Rxp)
    M[, u:(u + n_channels - 1)] <- do_iter(amp_matrix, k, N)
  }

  epsil <- 1 / sqrt(N) / 100

  if (epsil > tol) {
    message("Setting tolerance to ", round(epsil, 4))
    tol <- epsil
  }

  dim(M) <- c(n_channels, n_channels, n_lags)

  V <- JADE::rjd(M,
                 eps = tol,
                 maxiter = maxiter)$V

  ## create mixing matrix for output
  mixing_matrix <- MASS::ginv(Q, tol = 0) %*% V

  unmixing_matrix <- MASS::ginv(mixing_matrix, tol = 0)
  # rescale vecs
  scaling <- sqrt(colMeans(mixing_matrix^2))

  unmixing_matrix <- sweep(unmixing_matrix,
                           MARGIN = 1,
                           scaling,
                           `*`) # scaled weights

  mixing_matrix <- MASS::ginv(unmixing_matrix %*% diag(ncol(unmixing_matrix)),
                              tol = 0)

  dim(amp_matrix) <- c(n_channels,
                       n_times * n_epochs)

  S <- tcrossprod(unmixing_matrix,
                  as.matrix(data$signals))
  S <- as.data.frame(t(S))
  names(S) <- sprintf("Comp%03d", 1:ncol(S))


  list(M = mixing_matrix,
       W = unmixing_matrix,
       S = S)
}


#' Recreate channel timecourses from ICA decompositions.
#'
#' This function can be used to either recreate "mixed" (i.e. channel level)
#' timecourses from an ICA decomposition, or to apply a set of ICA weights to a
#' given dataset for the purpose of removing specific ICA components from that
#' dataset.
#'
#' @param data An \code{eeg_ICA} or \code{eeg_epochs} object.
#' @param ... Other parameters.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @export
apply_ica <- function(data, ...) {
  UseMethod("apply_ica", data)
}

#' @param comps Components to remove.
#' @describeIn apply_ica From given \code{eeg_ICA} object, recreate channel timecourses.
#' @export
apply_ica.eeg_ICA <- function(data,
                              comps = NULL,
                              ...) {
  ncomps <- ncol(data$mixing_matrix)
  new_mixmat <- data$mixing_matrix[1:(ncomps - 1)]
  new_mixmat[ , comps] <- 0
  new_dat <- mix(as.matrix(new_mixmat), as.matrix(data$signals))
  new_dat <- as.data.frame(t(new_dat))
  names(new_dat) <- data$chan_info$electrode
  out <-
    eeg_data(new_dat,
             srate = data$srate,
             events = data$events,
             chan_info = data$chan_info,
             timings = data$timings,
             reference = NULL,
             epochs = epochs(data))
  class(out) <- c("eeg_epochs", "eeg_data")
  out
}

#' @param decomp An \code{eeg_ICA} object.
#' @describeIn apply_ica Combine a specific set of ICA weights with any \code{eeg_epochs} object.
#' @export
apply_ica.eeg_epochs <- function(data,
                                 decomp,
                                 comps,
                                 ...) {
  if (!is.eeg_ICA(decomp)) {
    stop("An eeg_ICA object must be supplied as decomp.")
  }

  ncomps <- ncol(decomp$mixing_matrix) - 1
  nelecs <- ncol(data$signals)
  sources <- unmix(data$signals,
                   decomp$unmixing_matrix[, 1:nelecs])
  keep_comps <- names(decomp$mixing_matrix[, 1:ncomps])[-comps]
  new_sigs <- mix(sources[, -comps],
                  decomp$mixing_matrix[, keep_comps])
  colnames(new_sigs) <- data$chan_info$electrode
  data$signals <- tibble::as_tibble(new_sigs)
#
  data
}


mix <- function(sources,
                mixing_matrix) {
  t(tcrossprod(as.matrix(mixing_matrix), as.matrix(sources)))
}

unmix <- function(data,
                  unmixing_matrix) {
  t(tcrossprod(as.matrix(unmixing_matrix), as.matrix(data)))
}
