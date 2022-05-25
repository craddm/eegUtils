test_data <- import_raw("Newtest17-256.bdf")
test_data <- electrode_locations(test_data, montage = "biosemi64alpha",
                                 overwrite = TRUE)
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
                          n_freq = 10,
                          n_cycles = 3)
  expect_s3_class(tfr_test, "eeg_tfr")
})

test_that("Hanning method works", {
  tfr_test <- compute_tfr(demo_epochs,
                          method = "hanning",
                          foi = c(4, 30),
                          n_freq = 10,
                          n_cycles = 3)
  tfr_log <- compute_tfr(demo_epochs,
                         method = "hanning",
                         spacing = "log",
                         foi = c(4, 30),
                         n_freq = 10,
                         n_cycles = 3)
  expect_s3_class(tfr_test, "eeg_tfr")
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "hanning tfr plot",
    plot_tfr(tfr_test)
    )
  vdiffr::expect_doppelganger(
    "hanning tfr log-spaced plot",
    plot_tfr(tfr_log)
  )
})
