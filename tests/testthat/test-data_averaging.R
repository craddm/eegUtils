test_that("eeg_average produces sensible output", {
  demo_spat_evo <- eeg_average(demo_spatial)
  expect_equal(demo_spat_evo,
               readRDS("reference_files/demo_spatial_evoked.rds"))
  expect_s3_class(eeg_average(demo_spat_evo,
                              cols = "everything"),
                  "eeg_evoked")
})
