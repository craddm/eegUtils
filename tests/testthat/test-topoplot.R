context("topoplot functions")
load("EEGdat.rda")
demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)
demo_SOBI <- run_ICA(demo_epochs, pca = 10)

test_that("Topoplots", {
  skip_on_appveyor()
  skip_on_travis()
  vdiffr::expect_doppelganger("topoplot of epochs",
                              topoplot(demo_epochs,
                                       limits = c(-5, 5))
                              )
  vdiffr::expect_doppelganger("topoplot of epochs 150-200ms",
                              topoplot(demo_epochs,
                                       time_lim = c(.150, .200),
                                       limits = c(-5, 5))
                              )
  vdiffr::expect_doppelganger("GAM topo",
                              topoplot(EEGdat,
                                       time_lim = c(150, 200),
                                       method = "gam",
                                       limits = c(-3, 3))
                              )
  vdiffr::expect_doppelganger("ica topo",
                              topoplot(demo_SOBI,
                                       component = 1,
                                       time_lim = c(.12, .15),
                                       limits = c(-5, 5))
                              )
  vdiffr::expect_doppelganger("topoplot head limit",
                              topoplot(demo_epochs,
                                       time_lim = c(.150, .200),
                                       limits = c(-5, 5),
                                       interp_limit = "head")
                              )
})


