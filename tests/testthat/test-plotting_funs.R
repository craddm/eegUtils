
context("Plotting functions")
load("EEG_epochs.rda")
load("EEGdat.rda")
test_data <- import_raw("Newtest17-256.bdf")
demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

test_that("Plotting of single epoch timecourses works as expected", {
  vdiffr::expect_doppelganger("timecourse over all elecs",
                              plot_timecourse(EEGdat))
  vdiffr::expect_doppelganger("timecourse plot from EEGdat",
                              plot_timecourse(EEGdat,
                                              electrode = "Pz"))
  vdiffr::expect_doppelganger("baseline corr Pz", plot_timecourse(EEGdat,
                                              baseline = c(-100,
                                                           0),
                                              electrode = "Pz"))
})

test_that("Plotting of data with multiple epochs works as expected", {
  vdiffr::expect_doppelganger("epochs plot",
                              plot_timecourse(EEG_epochs))
  vdiffr::expect_doppelganger("Pz only epochs",
                              plot_timecourse(EEGdat,
                                              electrode = "Pz"))
  vdiffr::expect_doppelganger("Oz baseline corr epochs",
                              plot_timecourse(EEGdat,
                                              baseline = c(-200, 0),
                                              electrode = "Oz"))
})

test_that("Plotting of butterfly plots from epochs", {
  skip('skip')
  vdiffr::expect_doppelganger("butterfly epochs",
                              plot_butterfly(EEG_epochs))
})

test_that("Topoplots", {
  vdiffr::expect_doppelganger("topoplot of epochs",
                              topoplot(demo_epochs))
  vdiffr::expect_doppelganger("topoplot of epochs 150-200ms",
                              topoplot(EEG_epochs,
                                       time_lim = c(150, 200)))
  vdiffr::expect_doppelganger("GAM topo",
                              topoplot(EEGdat,
                                       time_lim = c(150, 200),
                                       method = "gam"))
})

test_that("erp_raster and erp_image function", {
  test_epo <- epoch_data(test_data, 255)
  expect_s3_class(erp_image(demo_epochs, electrode = "A13"), "gg")
  expect_s3_class(erp_raster(demo_epochs), "gg")

  expect_error(electrode_locations(test_epo, montage = "bio2"))
  test_epo <- electrode_locations(test_epo, montage = "biosemi64alpha")
  expect_is(test_epo$chan_info, "data.frame")

})
