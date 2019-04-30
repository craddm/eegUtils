#' Fit a GLM model to each timepoint for an individual subject.
#'
#' Fits a linear model to each timepoint using lm().
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param formula A regression formula for a GLM. See ?formula
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

#' @export
fit_glm.default <- function(formula,
                            data,
                            ...) {
  stop(paste("Objects of class", class(.data), "not currently supported"))

}

#' @export
fit_glm.eeg_epochs <- function(formula,
                               data,
                               family = NULL,
                               ...) {

  data <- tibble::as_tibble(as.data.frame(data,
                                          long = TRUE,
                                          coords = FALSE))
  data <- tidyr::nest(data,
                       -time,
                       -electrode,
                       .key = "signals"
                       )

  mdf <- stats::model.matrix(formula,
                             data$signals[[1]])
  lm_test <- function(mdf, x) {
    qr.coef(qr(mdf), x)
  }

  if (is.null(family)) {
    data$fit <-
    lapply(data$signals,
           function(x) stats::lm.fit(mdf,
                                     x$amplitude)$coefficients)
  } else {
    data$fit <-
      lapply(data$signals,
             function(x) stats::glm.fit(mdf,
                                       x$amplitude, family = family)$coefficients)
  }


  data$signals <- NULL
  data
}
