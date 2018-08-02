context("Test artefact rejection methods")

demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

test_that("FASTER runs correctly.", {
  test_FASTER <- eeg_FASTER(demo_epochs)
  expect_known_output(test_FASTER, "faster_t.Rdata")
  #expect_known_hash(test_FASTER, "04de56b5ba")
})
