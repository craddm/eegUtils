context("Frequency analysis plotting")

test_that("plot_psd produces a ggplot", {
  psd_test <- plot_psd(demo_epochs)
  vdiffr::expect_doppelganger("power spectrum for demo epochs", psd_test)
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 16,
                          n_cycles = 3)
  tfr_plot <- plot_tfr(tfr_test)
  vdiffr::expect_doppelganger("TFA plot for demo epochs", tfr_plot)

})
