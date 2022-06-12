test_that("simple plotting of gfp", {
  vdiffr::expect_doppelganger("gfp for demo_spatial",
                              plot_gfp(demo_spatial))
  vdiffr::expect_doppelganger("gfp for demo_spatial with colours",
                              plot_gfp(demo_spatial,
                                       cols = "epoch_labels"))
  vdiffr::expect_doppelganger("gfp for demo_spatial with trials kept",
                              plot_gfp(demo_spatial,
                                       cols = "epoch_labels",
                                       keep_trials = TRUE))
})
