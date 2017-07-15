#' Fit a glm model to each timepoint for an individual subject.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @param data An EEG dataset.
#' @importFrom purrr map
#' @import tidyr
#' @import dplyr
#' @import broom
#'


fit_glm_indiv <- function(data, dv, iv_1) {
  #options(contrasts=c('contr.sum','contr.poly'))
  if (!"epoch" %in% colnames(data)) {
    stop("Single trial data required for fitting.")
  }
  if (is_grouped_df(data)) {
    data <- ungroup(data)
  }
  data <- nest(data, condition, amplitude)
  data <- mutate(data, fit = map(data, ~tidy(lm(as.name(dv) ~ as.name(iv_1), data = .))))

}
