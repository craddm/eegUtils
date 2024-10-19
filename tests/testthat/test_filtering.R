test_data <- import_raw("Newtest17-256.bdf")

test_that("Filtering works for eeg_* objects", {

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

  test_epo <- epoch_data(test_data, 255)
  test_iir <- eeg_filter(test_epo, low_freq = 1, method = "iir")
  expect_s3_class(test_iir$signals, "tbl_df")
  expect_false(any(is.na(test_epo$signals[, 1])))

})


