#' Fit a glm model to each timepoint for an individual subject.
#'
#' @author Matt Craddock, \email{matt@@mattcraddock.com}
#' @param data An EEG dataset.
#' @param dv Column containing dependent variable. (e.g. amplitude)
#' @param iv_1 Column containing predictor.
#' @importFrom purrr map
#' @importFrom tidyr nest
#' @importFrom dplyr mutate ungroup is_grouped_df
#' @noRd

fit_glm_indiv <- function(data, dv, iv_1) {
  #options(contrasts=c('contr.sum','contr.poly'))
  if (!"epoch" %in% colnames(data)) {
    stop("Single trial data required for fitting.")
  }
  if (dplyr::is_grouped_df(data)) {
    data <- dplyr::ungroup(data)
  }
  data <- tidyr::nest(data, condition, amplitude)
  data <- dplyr::mutate(data,
                        fit = purrr::map(data,
                                         ~lm(as.name(dv) ~ as.name(iv_1),
                                                         data = .)))

}
