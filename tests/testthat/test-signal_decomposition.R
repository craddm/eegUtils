context("ssd test")
test_that("ssd runs", {
  skip_on_ci()
  ssd_demo <-
    eeg_decompose(demo_epochs,
                  sig_range = c(8, 12),
                  noise_range = c(7, 13))
  vdiffr::expect_doppelganger("topographical plot for ssd",
                              topoplot(ssd_demo,
                                       component = "Comp001",
                                       limits = c(-4, 3)))
})
