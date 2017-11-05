context("Plotting functions")

load("EEG_epochs.rda")
load("EEGdat.rda")


test_that("Plotting of single epoch timecourses works as expected", {
  expect_equal_to_reference(plot_timecourse(EEGdat), file = "avg_all_elecs.rds")
  expect_equal_to_reference(plot_timecourse(EEGdat, electrode = "Pz"),
                            file = "Pz_single_epoch.rds")
  expect_equal_to_reference(plot_timecourse(EEGdat, baseline = c(-100, 0),
                                            electrode = "Pz"),
                            file = "Pz_bl_corrected.rds")
})
