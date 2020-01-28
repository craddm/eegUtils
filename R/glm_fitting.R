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
#' @param formula An object of class "\code{formula}". Right-hand side A
#'   regression formula for a GLM. See ?formula and notes on use below
#' @param data An \code{eegUtils} object.
#' @param ... Any other arguments passed to (LM/GLM)
#' @importFrom tibble as_tibble
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
#'   as a baseline.
#' @describeIn fit_glm GLM fitting for \code{eeg_epochs}
#' @export
fit_glm.eeg_epochs <- function(formula,
                               data,
                               time_lim = NULL,
                               ...) {

  n_chans <- ncol(data$signals)
  n_epochs <- nrow(epochs(data))
  chan_names <- channel_names(data)

  # split the data timepoint-by-timepoint
  split_data <- split(
    data$signals,
    data$timings$time
  )

  n_times <- length(split_data)

  mod_terms <- all.vars(formula)

  if ("baseline" %in% mod_terms) {

    if (is.null(time_lim)) {
      stop("time_lim must be specified if using baseline fitting")
    }

    base_times <- get_epoch_baselines(
      data,
      time_lim
    )

    # basic design matrix with dummy baseline column
    design_mat <- stats::model.matrix(
      formula,
      cbind(epochs(data),
        baseline = 0
      )
    )

    # create one design matrix for each channel with the appropriate baseline
    # for each epoch
    mod_mats <- lapply(
      seq(1, n_chans),
      function(x) {
        model.matrix(
          formula,
          cbind(epochs(data),
            baseline = base_times[, x]
          )
        )
      }
    )

    n_preds <- ncol(mod_mats[[1]])

    # invert design matrices using qr/cholesky
    inverted_dm <- lapply(
      mod_mats,
      function(x) chol2inv(qr(x)$qr[1:n_preds, 1:n_preds])
    )

    # initialize empty lists for storing
    out <- vector(
      mode = "list",
      length = n_times
    )

    out_se <- vector(
      mode = "list",
      length = n_times
    )
    out_res <- vector(
      mode = "list",
      length = n_times
    )

    tmp <- vector(
      mode = "list",
      length = n_chans
    )

    tmp_se <- vector(
      mode = "list",
      length = n_chans
    )
    tmp_res <- vector(
      mode = "list",
      length = n_chans
    )

    resid_df <- n_epochs - n_preds

    tot_vars <- lapply(
      split_data,
      function(x) {
        colSums(scale(x,
          scale = FALSE
        )^2)
      }
    )

    # Outer loop through timepoints, inner loop through channels
    # Try other way to see if any faster
    for (timepoint in (seq_along(split_data))) {
      for (i in (seq_along(mod_mats))) {
        tmp_mod <- .lm.fit(
          mod_mats[[i]],
          split_data[[timepoint]][[i]]
        )
        tmp[[i]] <- coef(tmp_mod)
        tmp_res[[i]] <- sum(resid(tmp_mod)^2)
        tmp_se[[i]] <- sqrt(diag(inverted_dm[[i]] * (tmp_res[[i]] / resid_df)))
      }

      out[[timepoint]] <- do.call(
        cbind,
        tmp
      )
      out_se[[timepoint]] <- do.call(
        cbind,
        tmp_se
      )
      out_res[[timepoint]] <- do.call(
        cbind,
        tmp_res
      )
    }

    r_sq <- tibble::as_tibble(1 - do.call(
      rbind,
      out_res
    ) / do.call(
      rbind,
      tot_vars
    ))

    r_sq$time <- unique(data$timings$time)

  } else {

    design_mat <- stats::model.matrix(
      formula,
      epochs(data)
    )

    n_preds <- ncol(design_mat)
    # invert matrix using qr decomposition
    inverted_dm <- chol2inv(qr(design_mat)$qr[1:n_preds, 1:n_preds])
    resid_df <- nrow(design_mat) - n_preds

    tot_vars <- lapply(
      split_data,
      function(x) {
        colSums(scale(x,
          scale = FALSE
        )^2)
      }
    )

    out <- lapply(
      split_data,
      function(x) {
        get_lm(
          design_mat,
          as.matrix(x),
          inverted_dm,
          resid_df
        )
      }
    )

    out_se <- lapply(
      out,
      function(x) x$ses
    )

    res_vars <- lapply(
      out,
      function(x) x$res_vars
    )

    r_sq <- tibble::as_tibble(1 - do.call(
      rbind,
      res_vars
    ) / do.call(
      rbind,
      tot_vars
    ))

    out <- lapply(
      out,
      function(x) x$coefs
    )

    r_sq$time <- unique(data$timings$time)
  }

  all_coefs <- do.call(rbind, out)
  coef_names <- rep(
    colnames(design_mat),
    length(split_data)
  )

  all_coefs <- tibble::as_tibble(all_coefs)
  names(all_coefs) <- chan_names# channel_names(data)

  # for compatibility with other eegUtils structures, create a timings structure.
  timings <- tibble::tibble(
    time = rep(unique(data$timings$time),
      each = n_preds
    ),
    epoch = rep(
      1:n_preds,
      nrow(all_coefs) / n_preds
    )
  )

  epochs <-
    tibble::new_tibble(
      list(
        epoch = unique(timings$epoch),
        participant_id = rep(
          unique(epochs(data)$participant_id),
          n_preds
        ),
        recording = rep(
          unique(epochs(data)$participant_id),
          n_preds
        ),
        coefficient = unique(coef_names)
      ),
      nrow = n_preds,
      class = "epoch_info"
    )
  # all_coefs$coef <- coef_names #rep(coef_names, length(split_data))
  std_errs <- tibble::as_tibble(
    do.call(
      rbind,
      out_se
    )
  )

  names(std_errs) <- chan_names#channel_names(data)
  # std_errs$time <- rep(
  #   unique(data$timings$time),
  #   each = ncol(mdf)
  # )
  #
  # std_errs$coef <- coef_names # rep(coef_names, length(split_data))
  t_stats <- std_errs
  t_stats[, 1:n_chans] <- all_coefs[, 1:n_chans] / std_errs[, 1:n_chans]

  new_eeg_lm(
    "coefficients" = all_coefs,
    "std_err" = std_errs,
    "t_stats" = t_stats,
    "chan_info" = channels(data),
    "epochs" = epochs,
    "r_sq" = r_sq,
    "timings" = timings,
    "formula" = formula
  )
}

#' Calculate coefficients and standard errors
#'
#' @param mdf model matrix
#' @param x data to be modelled (matrix of trials x channels)
#' @param inverted_dm inverted matrix
#' @keywords internal
get_lm <- function(mdf,
                   x,
                   inverted_dm,
                   resid_df) {

  tmp_mod <- .lm.fit(mdf,
                     x)

  res_vars <- colSums(resid(tmp_mod)^2)

  out_se <- lapply(res_vars,
                   function(x) sqrt(diag(inverted_dm * (x/resid_df))))
  list("coefs" = coef(tmp_mod),
       "ses" = do.call(cbind,
                       out_se),
       "res_vars" = res_vars)
}

# get hat values for calculating robust standard errors
get_hat <- function(mdf,
                    inverted_dm) {
  diag(mdf %*% inverted_dm %*% t(mdf))
}

# robust standard error calculation
rob_se <- function(mdf,
                   u2,
                   inverted_dm) {
  ouchman <- sweep(mdf, 1, u2, "*")
  inverted_dm %*% (t(ouchman) %*% mdf) %*% inverted_dm
}


#
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
