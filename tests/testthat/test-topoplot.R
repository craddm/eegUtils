load("EEGdat.rda")
demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)
demo_SOBI <- run_ICA(demo_epochs, pca = 10)
demo_tfr <- compute_tfr(demo_epochs,
                        foi = c(4,30),
                        n_freq = 10)

test_that("Topoplots", {
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "topoplot of epochs",
    topoplot(demo_epochs,
             limits = c(-2.87, 4.69)))
  vdiffr::expect_doppelganger(
    "topoplot of epochs 150-200ms",
    topoplot(demo_epochs,
             time_lim = c(.150, .200),
             limits = c(-4, 4)))
  vdiffr::expect_doppelganger(
    "GAM topo",
    topoplot(EEGdat,
             time_lim = c(150, 200),
             method = "gam",
             limits = c(-2.25, 2.75)))
  vdiffr::expect_doppelganger(
    "head limit",
    topoplot(demo_epochs,
             time_lim = c(.15, .20),
             interp_limit = "head",
             limits = c(-3, 3)))
  vdiffr::expect_doppelganger(
    "multiple times",
    topoplot(demo_spatial,
             time_lim = list(0, .1, .2)))
})
