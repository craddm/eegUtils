context("Frequency analysis plotting")

test_that("plot_psd produces a ggplot", {
  psd_test <- plot_psd(demo_epochs)
  vdiffr::expect_doppelganger("power spectrum for demo epochs", psd_test)
})
