#' Fit a glm model to each timepoint for an individual subject.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data An EEG dataset.
#' @param formula A regression formula for a GLM. See ?formula
#' @importFrom purrr map
#' @importFrom tidyr nest
#' @importFrom dplyr mutate ungroup is_grouped_df
#' @noRd

fit_glm_indiv <- function(data, formula) {
  #options(contrasts=c('contr.sum','contr.poly'))
  if (!"epoch" %in% colnames(data)) {
    stop("Single trial data required for fitting.")
  }
  if (dplyr::is_grouped_df(data)) {
    data <- dplyr::ungroup(data)
  }
  data <- tidyr::nest(data,
                      condition,
                      amplitude)
  data <- dplyr::mutate(data,
                        fit = purrr::map(data,
                                         ~lm(formula,
                                             data = .)))

}
