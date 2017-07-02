#' Fit a glm model to each timepoint for an individual subject
#'
#' @param df
#'


fit_glm_indiv <- function(df, dv, iv_1) {
  options(contrasts=c('contr.sum','contr.poly'))

  split_factors <- list(time, electrode)

  df %>% nest(-time,-electrode) %>% mutate(fit = map(data, ~broom::tidy(lm(get(dv) ~ get(iv_1), data = .))))

}
