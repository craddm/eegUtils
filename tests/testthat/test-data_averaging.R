test_that("eeg_average produces sensible output", {
  demo_spat_evo <- eeg_average(demo_spatial)
  expect_s3_class(demo_spat_evo, "eeg_evoked")
})
