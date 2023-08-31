context("Test glm.R")

demo_tagged <-
  eegUtils::tag_epochs(
    eegUtils::tag_events(
      eegUtils::demo_epochs,
      c(208, 213, 215,
        207, 222, 219),
      event_label = c("Match", "Match", "Match",
                      "Mismatch", "Mismatch", "Mismatch")),
    event_label = c("Match",
                    "Mismatch"))

test_that("glm fitting works", {
  test_glm <- fit_glm(~event_label,
                      data = demo_tagged)
  expect_known_output(test_glm,
                      "test_glm.Rdata")
  test_glm_bl <- fit_glm(~event_label + baseline,
                         data = demo_tagged,
                         baseline = c(-.1, 0))
  expect_known_output(test_glm_bl,
                      "test_glm_bl.Rdata")
})
