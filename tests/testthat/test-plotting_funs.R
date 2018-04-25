
context("Plotting functions")
load("EEG_epochs.rda")
load("EEGdat.rda")
test_data <- import_raw("Newtest17-256.bdf")

test_that("Plotting of single epoch timecourses works as expected", {
  expect_equal_to_reference(plot_timecourse(EEGdat), file = "avg_all_elecs.rds")
  expect_equal_to_reference(plot_timecourse(EEGdat, electrode = "Pz"),
                            file = "Pz_single_epoch.rds")
  expect_equal_to_reference(plot_timecourse(EEGdat, baseline = c(-100, 0),
                                            electrode = "Pz"),
                            file = "Pz_bl_corrected.rds")
})

test_that("Plotting of data with multiple epochs works as expected", {
  expect_equal_to_reference(plot_timecourse(EEG_epochs), file = "avg_epochs_elecs.rds")
  expect_equal_to_reference(plot_timecourse(EEGdat, electrode = "Pz"),
                            file = "Pz_epochs.rds")
  expect_equal_to_reference(plot_timecourse(EEGdat, baseline = c(-200, 0),
                                            electrode = "Oz"),
                            file = "Pz_bl_epochs.rds")
})

test_that("Plotting of butterfly plots from epochs", {
  expect_equal_to_reference(plot_butterfly(EEG_epochs), file = "butterfly_epochs.rds")
})

test_that("Topoplots", {
  expect_equal_to_reference(topoplot(EEG_epochs), file = "topo_epochs.rds")
  expect_equal_to_reference(topoplot(EEG_epochs, time_lim = c(150, 200)),
                            file = "150_250_epochs.rds")
  expect_equal_to_reference(topoplot(EEGdat, time_lim = c(150, 200),
                                            method = "gam"),
                            file = "gam_topo.rds")
})

test_that("erp_raster and erp_image function", {
  test_epo <- epoch_data(test_data, 255)
  expect_s3_class(erp_image(test_epo, electrode = "A15"), "gg")
  expect_s3_class(erp_raster(test_epo), "gg")
})
