context("Plotting functions")
load("EEGdat.rda")
test_data <- import_raw("Newtest17-256.bdf")
demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

test_that("Plotting of data with multiple epochs works as expected", {
  vdiffr::expect_doppelganger("epochs plot",
                              plot_timecourse(demo_epochs))
  vdiffr::expect_doppelganger("A29 only epochs",
                              plot_timecourse(demo_epochs,
                                              electrode = "A29"))
  vdiffr::expect_doppelganger("A29 baseline corr epochs",
                              plot_timecourse(demo_epochs,
                                              baseline = c(-.2, 0),
                                              electrode = "A29"))
})

test_that("Plotting of butterfly plots from epochs", {
  vdiffr::expect_doppelganger("butterfly epochs",
                              plot_butterfly(demo_epochs))
  vdiffr::expect_doppelganger("butterfly epochs baseline",
                              plot_butterfly(demo_epochs,
                                             baseline = c(-.2, 0),
                                             electrode = "A29"))
})

test_that("Topoplots", {
  skip_on_appveyor()
  vdiffr::expect_doppelganger("topoplot of epochs",
                              topoplot(demo_epochs))
  vdiffr::expect_doppelganger("topoplot of epochs 150-200ms",
                              topoplot(demo_epochs,
                                       time_lim = c(.150, .200)))
  vdiffr::expect_doppelganger("GAM topo",
                              topoplot(EEGdat,
                                       time_lim = c(150, 200),
                                       method = "gam"))
})


test_that("erp_raster and erp_image function", {
  test_epo <- epoch_data(test_data, 255)
  expect_s3_class(erp_image(demo_epochs,
                            electrode = "A13"),
                  "gg")
  expect_s3_class(erp_raster(demo_epochs),
                  "gg")

  expect_error(electrode_locations(test_epo,
                                   montage = "bio2"))
  test_epo <- electrode_locations(test_epo,
                                  montage = "biosemi64alpha")
  expect_is(test_epo$chan_info,
            "data.frame")

})
