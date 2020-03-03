context("Test glm.R")

demo_tagged <- tag_epochs(tag_events(demo_epochs,
                                     c(208,
                                       213,
                                       215,
                                       207,
                                       222,
                                       219),
                                     event_label = c("Match",
                                                     "Match",
                                                     "Match",
                                                     "Mismatch",
                                                     "Mismatch",
                                                     "Mismatch")),
                          event_label = c("Match",
                                          "Mismatch"))

test_that("glm fitting works", {
  test_glm <- fit_glm(~event_label, data = demo_tagged)
  expect_known_output(test_glm,
                      "test_glm.Rdata")
})

