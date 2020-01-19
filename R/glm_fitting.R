#' Fit a linear model to EEG data
#'
#' Fits a linear model to each timepoint separately for each electrode.
#'
#' @section Notes On Usage:
#'
#'   The \code{fit_glm} function will fit a linear model to each timepoint for
#'   each electrode in the input dataset.
#'
#'   The formula is a standard R formula. Specify only the right-hand side of
#'   the formula i.e. any predictors.
#'
#'   The function allows flexible fitting of a baseline covariate, recognising
#'   the term \code{baseline} in the specified formula.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param formula An object of class "\code{formula}". Right-hand side A regression formula for a GLM. See ?formula
#'   and notes on use below
#' @param data An \code{eegUtils} object.
#' @param ... Any other arguments passed to (LM/GLM)
#' @importFrom purrr map
#' @importFrom tidyr nest
#' @importFrom dplyr mutate
#' @export

fit_glm <- function(formula,
                    data,
                    ...) {
  UseMethod("fit_glm", data)
}

fit_glm.default <- function(formula,
                            data,
                            ...) {
  stop(paste("Objects of class", class(.data), "not currently supported"))

}

#' @param time_lim Numeric vector of length 2 specifying time period to be used
#'   as a baseline
#' @describeIn fit_glm GLM fitting for \code{eeg_epochs}
#' @export
fit_glm.eeg_epochs <- function(formula,
                               data,
                               family = NULL,
                               time_lim = NULL,
                               ...) {

  zip <- split(data$signals,
               data$timings$time)

  mod_terms <- all.vars(formula)

  n_chans <- ncol(data$signals)

  if ("baseline" %in% mod_terms) {
    if (is.null(time_lim)) {
      stop("time_lim must be specified if using baseline fitting")
    }
    n_epochs <- nrow(epochs(data))

    base_times <- select_times(data,
                               time_lim = time_lim)
    base_times$signals <- as.matrix(base_times$signals)
    n_bl_times <- length(unique(base_times$timings$time))
    dim(base_times$signals) <- c(n_bl_times, n_epochs, n_chans)
    base_times <- colMeans(base_times$signals)

    #basic design matrix with dummy baseline column
    mdf <- stats::model.matrix(formula,
                               mutate(epochs(data),
                                      baseline = 0))
    n_preds <- ncol(mdf)
    # create one design matrix for each channel with the appropriate baseline
    # for each epoch
    mod_mats <- lapply(seq(1, n_chans),
                       function(x) model.matrix(formula,
                                                mutate(epochs(data),
                                                       baseline = base_times[, x])))
    std_errs <- lapply(mod_mats,
                       function(x) chol2inv(qr(x)$qr[1:n_preds, 1:n_preds]))
    #chol2inv(qr(mdf2)$qr[1:3, 1:3])
    out <- vector(mode = "list",
                  length = length(zip))

    out_se <- vector(mode = "list",
                     length = length(zip))
    out_res <- vector(mode = "list",
                     length = length(zip))
    tmp <- vector(mode = "list",
                  length = n_chans)
    tmp_se <- vector(mode = "list",
                     length = n_chans)
    tmp_res <- vector(mode = "list",
                      length = n_chans)
    n <- dim(mod_mats[[1]])[1]
    k <- dim(mod_mats[[1]])[2]
    resid_df <- n - k

    tot_vars <- lapply(zip,
                       function(x) colSums(scale(x,
                                                 scale = FALSE)^2))
    # TODO add r_sq with baseline

    for (n_epo in (seq(1, length(zip)))) {
      for (i in (seq(1, n_chans))) {
        tmp_mod <- .lm.fit(mod_mats[[i]],
                           zip[[n_epo]][[i]])
        tmp[[i]] <- coef(tmp_mod)
        tmp_res[[i]] <- sum(resid(tmp_mod)^2)
        tmp_se[[i]] <- sqrt(diag(std_errs[[i]] * (tmp_res[[i]]/resid_df)))
      }
      out[[n_epo]] <- do.call(cbind,
                              tmp)
      out_se[[n_epo]] <- do.call(cbind,
                                 tmp_se)
      out_res[[n_epo]] <- do.call(cbind,
                                  tmp_res)
    }

    r_sq <- tibble::as_tibble(1 - do.call(rbind,
                                          out_res) / do.call(rbind,
                                                             tot_vars))
    r_sq$time <- unique(data$timings$time)

  } else {

    mdf <- stats::model.matrix(formula,
                               epochs(data))
    #std_errs <- solve(crossprod(mdf))
    n_preds <- ncol(mdf)
    #invert matrix using qr decomposition
    std_errs <- chol2inv(qr(mdf)$qr[1:n_preds, 1:n_preds])
    resid_df <- nrow(mdf) - ncol(mdf)
    tot_vars <- lapply(zip,
                       function(x) colSums(scale(x,
                                                 scale = FALSE)^2))

    # now_t <- function(mdf, x, std_errs) {
    #   tmp_mod <- .lm.fit(mdf, x)
    #
    #   # ses <-  apply(resid(tmp_mod),
    #   #               2,
    #   #               crossprod) #/ resid_df
    #   res_vars <- resid(tmp_mod)^2
    #
    #   XDXs <- apply(res_vars,
    #                 2,
    #                 function(x) diag(rob_sd(mdf, x, std_errs)))
    #
    #   out_se <- sqrt(XDXs)
    #   list("coefs" = coef(tmp_mod),
    #        "ses" = out_se,
    #        "res_vars" = colSums(res_vars))
    # }
    out <- lapply(zip,
                  function(x) new_t(mdf,
                                    as.matrix(x),
                                    std_errs,
                                    resid_df))
    out_se <- lapply(out,
                     function(x) x$ses)
    res_vars <- lapply(out,
                       function(x) x$res_vars)
    r_sq <- tibble::as_tibble(1 - do.call(rbind,
                                          res_vars) / do.call(rbind,
                                                              tot_vars))
    out <- lapply(out,
                  function(x) x$coefs)
    r_sq$time <- unique(data$timings$time)
  }

  all_coefs <- do.call(rbind, out)
  coef_names <- rep(colnames(mdf),
                    length(zip))
  all_coefs <- tibble::as_tibble(all_coefs)
  names(all_coefs) <- channel_names(data)
  all_coefs$time <- rep(unique(data$timings$time),
                        each = ncol(mdf))
  all_coefs$coef <- coef_names #rep(coef_names, length(zip))
  std_errs <- tibble::as_tibble(do.call(rbind,
                                        out_se))
  names(std_errs) <- channel_names(data)
  std_errs$time <- rep(unique(data$timings$time),
                       each = ncol(mdf))
  std_errs$coef <- coef_names #rep(coef_names, length(zip))
  t_stats <- std_errs
  t_stats[, 1:n_chans] <- all_coefs[, 1:n_chans] / std_errs[, 1:n_chans]

  list("coefficients" = all_coefs,
       "std_err" = std_errs,
       "t_stats" = t_stats,
       "chan_info" = channels(data),
       "r_sq" = r_sq)
}

#' @param mdf model matrix
#' @param x data to be modelled (matrix of trials x channels)
#' @param std_errs
new_t <- function(mdf,
                  x,
                  std_errs,
                  resid_df) {

  tmp_mod <- .lm.fit(mdf,
                     x)

  res_vars <- colSums(resid(tmp_mod)^2)

  out_se <- lapply(res_vars,
                   function(x) sqrt(diag(std_errs * (x/resid_df))))
  list("coefs" = coef(tmp_mod),
       "ses" = do.call(cbind,
                       out_se),
       "res_vars" = res_vars)
}

get_hat <- function(mdf, std_errs) {
  diag(mdf %*% std_errs %*% t(mdf))
}

rob_se <- function(mdf, u2, XX) {
  ouchman <- sweep(mdf, 1, u2, "*")
  XX %*% (t(ouchman) %*% mdf) %*% XX
}
