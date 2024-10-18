#' Independent Component Analysis for EEG data
#'
#' Performs Independent Component Analysis for electroencephalographic data.
#' Currently only available with on epoched data. Implements three different
#' methods of ICA - 'fastica', 'extended Infomax', and 'Second-Order Blind
#' Identification (SOBI)'. The resulting `eeg_ICA` objects can be used largely
#' like `eeg_epochs` objects.
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
#'   In addition, PCA will be required if the data is not full rank. This is
#'   typical when using average reference, when the data rank will be
#'   n_electrodes - 1.
#'
#' @param data Data to be ICAed.
#' @param ... Other parameters passed to function.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @return An `eeg_ICA` object containing an ICA decomposition
#' @family decompositions
#' @examples
#' sobi_demo <-
#'   run_ICA(demo_epochs,
#'           pca = 10)
#'  sobi_demo
#'  # We can plot the resulting spatial filters using `topoplot()`
#'  topoplot(sobi_demo, 1:2)
#'  \dontrun{ view_ica(sobi_demo) }
#' @export

run_ICA <- function(data, ...) {
  UseMethod("run_ICA", data)
}

#' @export
run_ICA.default <- function(data,
                            ...) {
  stop("Not supported for objects of this class.")
}

#' @param method "sobi" (default), "fastica", "infomax", or "imax". "infomax"
#'   uses the implementation from the `ica` package, whereas `imax` uses the
#'   implementation from the `infomax` package, which is based on the `EEGLAB`
#'   implementation.
#' @param maxit Maximum number of iterations of the Infomax and Fastica ICA
#'   algorithms.
#' @param tol Convergence tolerance for fastica and infomax. Defaults to 1e-06.
#' @param pca Reduce the number of dimensions using PCA before running ICA.
#'   Numeric,  >1 and < number of channels
#' @param centre Defaults to TRUE. Centre the data on zero by subtracting the
#'   column mean. See notes on usage.
#' @param alg Use "gradient descent" or "newton" algorithm for extended infomax.
#'   Defaults to "gradient". Ignored if method != "infomax".
#' @param rateanneal Annealing rate for extended infomax. Ignored if method !=
#'   "infomax".
#' @param rate Learning rate for extended infomax. Ignored if method !=
#'   "infomax".
#' @param verbose Print informative messages to console.
#' @param return "full" or "weights". "full" returns the mixing and unmixing
#'   matrices and the source timecourses. "weights" returns only the mixing and
#'   unmixing matrices. Defaults to "full".
#' @param options List of option
#' @describeIn run_ICA Run ICA on an `eeg_epochs` object
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
                               verbose = TRUE,
                               return = c("full", "weights"),
                               params = NULL,
                               ...) {

  orig_chans <- channel_names(data)

  pca_result <- perform_pca(data, pca)
  data <- pca_result$data
  pca_flag <- pca_result$pca_flag
  pca_decomp <- pca_result$pca_decomp

  rank_check <- qr(data$signals)$rank

  if (rank_check < ncol(data$signals)) {
    stop(paste(
      "Data is rank deficient. PCA required. Detected rank: ",
      rank_check
    ))
  }

  if (!(method %in% c("sobi", "fastica", "infomax", "fica", "imax"))) {
    stop("Unknown method; available methods are sobi, fastica, infomax, fica, and imax.")
  }

  # Run ICA
  ICA_out <- run_ica_method(data, method, maxit, tol, centre, alg, rateanneal, rate, verbose)

  # Process ICA output
  processed_output <- process_ica_output(ICA_out, method, pca, pca_flag, pca_decomp, orig_chans)

  create_eeg_ica_object(processed_output, data, method, tol, maxit, return)
}

#' @export
run_ICA.eeg_group <- function(data,
                              ...) {
  stop("Not supported for group data.")
}

#' SOBI ICA
#'
#' Internal function for running SOBI ICA on an `eeg_epochs` object
#'
#' @param data Data to be ICAed.
#' @param maxiter Maximum number of iterations of the joint diagonalization
#' @param tol convergence tolerance.
#' @param pca Number of PCA components.
#' @param centre Mean centre signals.
#' @param verbose Print informative messages.
#' @param whitening Defaults to pca, options are pca or zca
#' @author A. Belouchrani and A. Cichocki. Adapted to R by Matt Craddock
#'   \email{matt@@mattcraddock.com}
#' @keywords internal

sobi_ICA <- function(data,
                     maxiter,
                     tol,
                     pca,
                     centre,
                     verbose = TRUE,
                     whitening = "pca") {

  n_epochs <- nrow(epochs(data))
  n_channels <- ncol(data$signals)
  n_times <- length(unique(data$timings$time))

  # Optimize: Use integer division for n_lags
  n_lags <- min(100, n_times %/% 3)

  if (centre) {
    data <- rm_baseline(data, verbose = FALSE)
  }

  cov_matrix <- cov(as.matrix(data$signals))
  eigen_decomp <- eigen(cov_matrix)
  if (identical(whitening, "pca")) {
    Q <- diag(1 / sqrt(eigen_decomp$values)) %*% t(eigen_decomp$vectors)
  } else {
    Q <- eigen_decomp$vectors %*% diag(1 / sqrt(eigen_decomp$values)) %*% t(eigen_decomp$vectors)
  }
  # Reshape amp_matrix to have three dimensions
  amp_matrix <- array(Q %*% t(as.matrix(data$signals)), dim = c(n_channels, n_times, n_epochs))

  # Preallocate matrix M
  M <- array(0, dim = c(n_channels, n_channels, n_lags))

  # Use do_iter function for cross-covariance calculation
  for (k in 1:n_lags) {
    M[, , k] <- do_iter(amp_matrix, k - 1, n_times)
  }

  epsil <- 1 / sqrt(n_times) / 100
  tol <- max(tol, epsil)

  V <- JADE::rjd(M, eps = tol, maxiter = maxiter)$V

  mixing_matrix <- solve(Q) %*% V

  var_order <- order(colSums(mixing_matrix^2), decreasing = TRUE)
  mixing_matrix <- mixing_matrix[, var_order]

  unmixing_matrix <- solve(mixing_matrix)

  S <- t(unmixing_matrix %*% t(as.matrix(data$signals)))
  S <- as.data.frame(S)
  names(S) <- sprintf("Comp%03d", seq_len(ncol(S)))

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
#' @param data An `eeg_ICA` or `eeg_epochs` object.
#' @param ... Other parameters.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @examples
#' test_ica <- run_ICA(demo_epochs, pca = 10)
#' plot_butterfly(demo_epochs)
#' # Reconstruct the original data from the ICA decomposition.
#' # Note that the ICA process subtracts the mean from each epoch,
#' # so the reconstructed plot may look slightly different to the original.
#' plot_butterfly(apply_ica(test_ica))
#' # Remove component 2 from the data
#' plot_butterfly(apply_ica(demo_epochs, test_ica, comps = 2))
#' @export
apply_ica <- function(data, ...) {
  UseMethod("apply_ica", data)
}

#' @param comps Components to remove.
#' @describeIn apply_ica From given `eeg_ICA` object, recreate channel
#'   timecourses.
#' @export
apply_ica.eeg_ICA <- function(data,
                              comps = NULL,
                              ...) {
  ncomps <- ncol(data$mixing_matrix)
  new_mixmat <- data$mixing_matrix[1:(ncomps - 1)]
  new_mixmat[, comps] <- 0
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

#' @param decomp An `eeg_ICA` object.
#' @describeIn apply_ica Combine a specific set of ICA weights with any
#'   `eeg_epochs` object.
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

  if (is.character(comps)) {
    comps_missing <- which(!(comps %in% channel_names(decomp)))
    keep_comps <- which(!channel_names(decomp) %in% comps)
  } else {
    keep_comps <- seq(1, ncomps)[-comps]
  }

  sources <- unmix(data$signals,
                   decomp$unmixing_matrix[, 1:nelecs])
  new_sigs <- mix(sources[, keep_comps],
                  decomp$mixing_matrix[, keep_comps])
  colnames(new_sigs) <- data$chan_info$electrode
  data$signals <- tibble::as_tibble(new_sigs)
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

vaf_mix <- function(mixing_matrix) {
  comp_var <- colSums(mixing_matrix^2)
  comp_var / sum(comp_var)
}

get_vaf <- function(decomp) {
  ncomps <- ncol(decomp$mixing_matrix) - 1
  comp_var <- colSums(decomp$mixing_matrix[, 1:ncomps]^2)
  comp_var / sum(comp_var)
}


perform_pca <- function(data, pca) {
  if (!is.null(pca)) {
    stopifnot("`pca` must be numeric" = is.numeric(pca))
    pca <- floor(pca)
    message("Reducing data to ", pca, " dimensions using PCA.")
    pca_decomp <- eigen(stats::cov(data$signals))$vectors
    data$signals <- as.data.frame(as.matrix(data$signals) %*% pca_decomp[, 1:pca])
    pca_flag <- TRUE
  } else {
    pca <- ncol(data$signals)
    pca_flag <- FALSE
    pca_decomp <- NULL
  }
  list(data = data, pca_flag = pca_flag, pca_decomp = pca_decomp)
}

run_ica_method <- function(data, method, maxit, tol, centre, alg, rateanneal, rate, verbose) {
  switch(method,
    sobi = run_sobi_ica(data, maxit, tol, centre),
    fastica = run_fastica(data, maxit, tol, centre),
    infomax = run_infomax(data, maxit, tol, centre, alg, rateanneal, rate),
    fica = run_fica(data, maxit, tol),
    imax = run_imax(data, maxit, tol, centre, verbose),
    stop("Invalid method")
  )
}

run_sobi_ica <- function(data, maxit, tol, centre) {
  if (!requireNamespace("JADE", quietly = TRUE)) {
    stop("Package \"JADE\" needed to use SOBI. Please install it.", call. = FALSE)
  }
  message("Running SOBI ICA.")
  sobi_ICA(data, maxiter = maxit, tol = tol, pca = NULL, centre = centre)
}

run_fastica <- function(data, maxit, tol, centre) {
  if (!requireNamespace("ica", quietly = TRUE)) {
    stop("Package \"ica\" needed to use fastica. Please install it.", call. = FALSE)
  }
  message("Running fastica (ica).")
  ica::icafast(data$signals, nc = ncol(data$signals), maxit = maxit, tol = tol, center = centre)
}

run_infomax <- function(data, maxit, tol, centre, alg, rateanneal, rate) {
  if (!requireNamespace("ica", quietly = TRUE)) {
    stop("Package \"ica\" needed to use infomax. Please install it.", call. = FALSE)
  }
  message("Running extended-Infomax (ica).")
  ica::icaimax(data$signals,
    nc = ncol(data$signals), maxit = maxit, fun = "ext",
    tol = tol, center = centre, alg = alg, rateanneal = rateanneal, rate = rate
  )
}

run_fica <- function(data, maxit, tol) {
  if (!requireNamespace("fICA", quietly = TRUE)) {
    stop("Package \"fICA\" needed to use the fICA implementation of fastica. Please install it.", call. = FALSE)
  }
  message("Running fastica (fICA).")
  fICA::fICA(as.matrix(data$signals), eps = tol, maxiter = maxit, method = "sym2")
}

run_imax <- function(data, maxit, tol, centre, verbose) {
  if (!requireNamespace("infomax", quietly = TRUE)) {
    stop("Package \"infomax\" needed to use the \"imax\" method, please install it from https://github.com/eegverse/infomax")
  }
  message("Running extended-Infomax (infomax).")
  infomax::run_infomax(data$signals, tol = tol, centre = centre, maxiter = maxit, whiten = "sqrtm", verbose = verbose)
}

process_ica_output <- function(ICA_out, method, pca, pca_flag, pca_decomp, orig_chans) {
  if (method == "fica") {
    process_fica_output(ICA_out, pca, pca_flag, pca_decomp, orig_chans)
  } else {
    process_other_ica_output(ICA_out, pca, pca_flag, pca_decomp, orig_chans)
  }
}

process_fica_output <- function(ICA_out, pca, pca_flag, pca_decomp, orig_chans) {
  ICA_out$W <- ICA_out$W[, 1:pca]
  mixing_matrix <- MASS::ginv(ICA_out$W, tol = 0)

  if (pca_flag) {
    mixing_matrix <- pca_decomp[, 1:pca] %*% mixing_matrix
  }

  var_order <- order(vaf_mix(mixing_matrix), decreasing = TRUE)
  mixing_matrix <- mixing_matrix[, var_order]

  unmixing_matrix <- as.data.frame(MASS::ginv(mixing_matrix, tol = 0))
  mixing_matrix <- as.data.frame(mixing_matrix)

  names(mixing_matrix) <- sprintf("Comp%03d", 1:pca)
  mixing_matrix$electrode <- orig_chans
  names(unmixing_matrix) <- orig_chans
  unmixing_matrix$Component <- sprintf("Comp%03d", 1:pca)

  ICA_out$S <- ICA_out$S[, var_order]
  colnames(ICA_out$S) <- sprintf("Comp%03d", 1:pca)
  ICA_out$S <- tibble::as_tibble(ICA_out$S[, 1:pca])

  list(mixing_matrix = mixing_matrix, unmixing_matrix = unmixing_matrix, S = ICA_out$S)
}

process_other_ica_output <- function(ICA_out, pca, pca_flag, pca_decomp, orig_chans) {
  ICA_out$S <- as.data.frame(ICA_out$S)
  names(ICA_out$S) <- sprintf("Comp%03d", seq_len(ncol(ICA_out$S)))

  if (pca_flag) {
    mixing_matrix <- as.data.frame(pca_decomp[, 1:pca] %*% ICA_out$M)
    unmixing_matrix <- as.data.frame(MASS::ginv(as.matrix(mixing_matrix), tol = 0))
    names(mixing_matrix) <- sprintf("Comp%03d", 1:pca)
    names(unmixing_matrix) <- orig_chans
    mixing_matrix$electrode <- orig_chans
    unmixing_matrix$Component <- sprintf("Comp%03d", 1:pca)
  } else {
    mixing_matrix <- as.data.frame(ICA_out$M)
    names(mixing_matrix) <- sprintf("Comp%03d", seq_len(ncol(ICA_out$M)))
    mixing_matrix$electrode <- orig_chans
    unmixing_matrix <- as.data.frame(ICA_out$W)
    names(unmixing_matrix) <- orig_chans
    unmixing_matrix$Component <- sprintf("Comp%03d", seq_len(ncol(ICA_out$S)))
  }

  list(mixing_matrix = mixing_matrix, unmixing_matrix = unmixing_matrix, S = ICA_out$S)
}

create_eeg_ica_object <- function(processed_output, data, method, tol, maxit, return) {
  eeg_ICA(
    mixing_matrix = processed_output$mixing_matrix,
    unmixing_matrix = processed_output$unmixing_matrix,
    signals = tibble::as_tibble(processed_output$S),
    timings = data$timings,
    events = data$events,
    chan_info = data$chan_info,
    srate = data$srate,
    epochs = data$epochs,
    algorithm = list(
      algorithm = method,
      tol = tol,
      max_iterations = maxit
    ),
    contents = return
  )
}
