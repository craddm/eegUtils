test_data <- import_raw("Newtest17-256.bdf")

test_that("IIR filtering works for eeg_data objects", {

  test_iir <- eeg_filter(test_data, low_freq = 1, method = "iir")
  expect_s3_class(test_iir, "eeg_data")
  expect_s3_class(test_iir$signals, "tbl_df")
  test_iir <- eeg_filter(test_data, low_freq = 1, high_freq = 40, method = "iir")
  expect_s3_class(test_iir, "eeg_data")
  expect_s3_class(test_iir$signals, "tbl_df")
  test_iir <- eeg_filter(test_data, high_freq = 40, method = "iir")
  expect_s3_class(test_iir, "eeg_data")
  expect_s3_class(test_iir$signals, "tbl_df")
  test_iir <- eeg_filter(test_data, low_freq = 40, high_freq = 30, method = "iir")
  expect_s3_class(test_iir, "eeg_data")
  expect_s3_class(test_iir$signals, "tbl_df")

})

test_that("FIR filtering works for eeg_data objects", {

  test_fir <- eeg_filter(test_data, low_freq = 1, method = "fir")
  expect_s3_class(test_fir, "eeg_data")
  expect_s3_class(test_fir$signals, "tbl_df")
  expect_equal(nrow(test_fir$signals), nrow(test_data$signals))

  test_fir <- eeg_filter(test_data, high_freq = 40, method = "fir")
  expect_s3_class(test_fir, "eeg_data")
  expect_s3_class(test_fir$signals, "tbl_df")

  test_fir <- eeg_filter(test_data, low_freq = 1, high_freq = 40, method = "fir")
  expect_s3_class(test_fir, "eeg_data")
  expect_s3_class(test_fir$signals, "tbl_df")

})

test_that("Epoch-aware IIR filtering works for eeg_epochs", {

  test_epo <- epoch_data(test_data, 255)
  test_iir <- eeg_filter(test_epo, low_freq = 1, method = "iir")
  expect_s3_class(test_iir, "eeg_epochs")
  expect_s3_class(test_iir$signals, "tbl_df")
  expect_false(any(is.na(test_iir$signals[, 1])))
  expect_equal(nrow(test_iir$signals), nrow(test_epo$signals))

})

test_that("Epoch-aware FIR filtering works for eeg_epochs", {

  test_epo <- epoch_data(test_data, 255)
  # FIR with band-pass should work and preserve epoch count
  test_fir <- eeg_filter(test_epo, low_freq = 1, high_freq = 40, method = "fir")
  expect_s3_class(test_fir, "eeg_epochs")
  expect_s3_class(test_fir$signals, "tbl_df")
  expect_equal(nrow(test_fir$signals), nrow(test_epo$signals))
  expect_false(any(is.na(test_fir$signals[, 1])))

})

test_that("Filtered signals have no NAs", {

  test_fir <- eeg_filter(test_data, high_freq = 40, method = "fir")
  expect_false(any(is.na(test_fir$signals)))

  test_iir <- eeg_filter(test_data, high_freq = 40, method = "iir")
  expect_false(any(is.na(test_iir$signals)))

})

test_that("Filtered signal length is preserved", {

  test_fir <- eeg_filter(test_data, high_freq = 40, method = "fir")
  expect_equal(nrow(test_fir$signals), nrow(test_data$signals))

  test_iir <- eeg_filter(test_data, high_freq = 40, method = "iir")
  expect_equal(nrow(test_iir$signals), nrow(test_data$signals))

})

test_that("FIR low-pass attenuates high frequencies", {

  # Filter at 30 Hz then check the PSD is very low above 40 Hz
  test_fir <- eeg_filter(test_data, high_freq = 30, method = "fir")
  psd_result <- compute_psd(test_fir, n_fft = 256, keep_all = FALSE)
  # Average power across all channels (exclude frequency column)
  chan_cols <- setdiff(names(psd_result), "frequency")
  avg_power <- rowMeans(psd_result[, chan_cols])
  high_freq_power <- mean(avg_power[psd_result$frequency > 50])
  low_freq_power <- mean(avg_power[psd_result$frequency > 5 &
                                     psd_result$frequency < 20])
  expect_true(high_freq_power < low_freq_power)

})
