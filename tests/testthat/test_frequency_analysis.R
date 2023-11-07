test_data <- import_raw("Newtest17-256.bdf")
test_data <- electrode_locations(test_data, montage = "biosemi64alpha",
                                 overwrite = TRUE)
test_epo <- epoch_data(test_data, 255)

test_that("PSD computation runs correctly.", {
  expect_snapshot(compute_psd(test_data))
  expect_snapshot(compute_psd(test_epo))
})

test_that("TFA works on epoched data with and without trimmed edges", {
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 10,
                          n_cycles = 3)
  expect_s3_class(tfr_test, "eeg_tfr")
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 10,
                          n_cycles = 3,
                          trim_edges = FALSE)
  tfr_trials <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 10,
                          n_cycles = 3,
                          keep_trials = TRUE)
  expect_s3_class(tfr_trials, "eeg_tfr")
  expect_identical(dim(tfr_trials$signals)[[1]],
                   nrow(epochs(tfr_trials)))
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "untrimmed edges",
    plot_tfr(tfr_test))
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
