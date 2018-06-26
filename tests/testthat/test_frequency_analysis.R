context("Test frequency analyses")

test_data <- import_raw("Newtest17-256.bdf")
test_data <- electrode_locations(test_data, montage = "biosemi64alpha")
test_epo <- epoch_data(test_data, 255)

test_that("PSD computation runs correctly.", {
  test_psd <- compute_psd(test_data)
  expect_known_hash(test_psd, "f6c17ba5ef")
  test_epo_psd <- compute_psd(test_epo)
  expect_equal(length(unique(test_epo_psd$epoch)), 39)
  expect_equal(length(unique(test_epo_psd$frequency)), 129)
})
