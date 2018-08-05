context("Test frequency analyses")

test_data <- import_raw("Newtest17-256.bdf")
test_data <- electrode_locations(test_data, montage = "biosemi64alpha")
test_epo <- epoch_data(test_data, 255)
tmp_psd <- "psd_output.Rdata"

test_that("PSD computation runs correctly.", {
  test_psd <- compute_psd(test_data)
  expect_known_output(test_psd, tmp_psd)
  test_epo_psd <- compute_psd(test_epo)
  expect_equal(length(unique(test_epo_psd$epoch)), 39)
  expect_equal(length(unique(test_epo_psd$frequency)), 128)
})

test_that("TFA works on epoched data", {
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 16,
                          n_cycles = 3)
  expect_s3_class(tfr_test, "eeg_tfr")
})
