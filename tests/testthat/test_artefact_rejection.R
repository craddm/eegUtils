context("Test artefact rejection methods")

test_data <- import_raw("Newtest17-256.bdf")
test_data <- electrode_locations(test_data, montage = "biosemi64alpha")
test_epo <- epoch_data(test_data, 255)

test_that("FASTER runs correctly.", {
  #test_FASTER <- eeg_FASTER(test_epo)
  #expect_known_hash(test_FASTER, "04de56b5ba")
})
