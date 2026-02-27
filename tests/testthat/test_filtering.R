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

# fir1 has a constant offset of ~-0.00097 in normalised frequency (W = f/Nyquist),
# so the -6 dB crossing falls slightly below the requested cutoff. In Hz the
# offset scales linearly with sample rate: ~-0.06 Hz at 128 Hz, ~-0.12 Hz at
# 256 Hz, ~-0.25 Hz at 512 Hz, ~-0.48 Hz at 1000 Hz. The offset is the same
# regardless of window type. Tests below use srate = 256, so +-0.5 Hz tolerance
# comfortably accommodates the ~-0.12 Hz offset seen there.

test_that("FIR low-pass -6dB crossing is within 0.5 Hz of specified cutoff", {
  resp <- filter_response(high_freq = 30, srate = 256)
  cutoff_idx <- which(resp$magnitude <= 0.5)[1]
  expect_true(!is.na(cutoff_idx))
  expect_true(abs(resp$frequency[cutoff_idx] - 30) < 0.5)
})

test_that("FIR high-pass -6dB crossing is within 0.5 Hz of specified cutoff", {
  resp <- filter_response(low_freq = 1, srate = 256)
  cutoff_idx <- which(resp$magnitude >= 0.5)[1]
  expect_true(!is.na(cutoff_idx))
  expect_true(abs(resp$frequency[cutoff_idx] - 1) < 0.5)
})

test_that("FIR bandpass magnitude at centre frequency is near 0dB", {
  resp <- filter_response(low_freq = 8, high_freq = 30, srate = 256)
  centre_freq <- sqrt(8 * 30)  # geometric mean ~ 15.5 Hz
  centre_idx <- which.min(abs(resp$frequency - centre_freq))
  expect_true(resp$magnitude_db[centre_idx] > -3)
})

test_that("IIR low-pass -3dB crossing is within 0.5 Hz of specified cutoff", {
  resp <- filter_response(high_freq = 30, srate = 256, method = "iir")
  cutoff_idx <- which(resp$magnitude <= (1 / sqrt(2)))[1]
  expect_true(!is.na(cutoff_idx))
  expect_true(abs(resp$frequency[cutoff_idx] - 30) < 0.5)
})

test_that("hann window filter_order auto works and cutoff is correct", {
  resp <- filter_response(high_freq = 30, srate = 256, window = "hann")
  cutoff_idx <- which(resp$magnitude <= 0.5)[1]
  expect_true(!is.na(cutoff_idx))
  expect_true(abs(resp$frequency[cutoff_idx] - 30) < 0.5)
})

test_that("blackman window filter_order auto works and cutoff is correct", {
  resp <- filter_response(high_freq = 30, srate = 256, window = "blackman")
  cutoff_idx <- which(resp$magnitude <= 0.5)[1]
  expect_true(!is.na(cutoff_idx))
  expect_true(abs(resp$frequency[cutoff_idx] - 30) < 0.5)
})
