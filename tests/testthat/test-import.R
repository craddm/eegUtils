context("Test import of bdf")

test_data <- import_raw("Newtest17-256.bdf")

test_that("Import and epoching of bdf files works correctly", {

  skip_on_cran()

  expect_s3_class(test_data, "eeg_data")
  expect_equal(test_data, import_raw("Newtest17-256.bdf", fast_bdf = FALSE))
  expect_true(is.eeg_data(test_data))
  expect_false(is.eeg_epochs(test_data))
  expect_true(is.data.frame(as.data.frame(test_data)))
  expect_true(is.data.frame(as.data.frame(test_data,
                                          long = TRUE)))
  expect_identical(names(test_data$signals),
                   c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
                     "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16"))

  test_epo <- epoch_data(test_data, 255)
  expect_s3_class(test_epo, "eeg_epochs")
  expect_true(all(c("epoch", "time") %in% names(test_epo$events)))
  expect_equal(length(unique(test_epo$timings$epoch)),
               length(unique(test_epo$events$epoch)))

  test_evo <- eeg_average(test_epo)
  expect_s3_class(test_evo, "eeg_evoked")

})



test_that("Filtering works for eeg_* objects", {

  skip_on_cran()

  test_iir <- iir_filt(test_data, 1)
  expect_s3_class(test_iir, "eeg_data")
  expect_is(test_iir$signals, "tbl_df")
  test_iir <- iir_filt(test_data, 1, 40)
  expect_s3_class(test_iir, "eeg_data")
  expect_is(test_iir$signals, "tbl_df")
  test_iir <- iir_filt(test_data, high_freq = 40)
  expect_s3_class(test_iir, "eeg_data")
  expect_is(test_iir$signals, "tbl_df")
  test_iir <- iir_filt(test_data, low_freq = 40, high_freq = 30)
  expect_s3_class(test_iir, "eeg_data")
  expect_is(test_iir$signals, "tbl_df")

  test_epo <- epoch_data(test_data, 255)
  test_iir <- iir_filt(test_epo, 1)
  expect_is(test_iir$signals, "tbl_df")
  expect_false(any(is.na(test_epo$signals[, 1])))

})

test_that("Interpolation works for eeg_* objects", {

  skip_on_cran()

  expect_error(interp_elecs(test_data, "A2"))
  test_elecs <- electrode_locations(test_data, montage = "biosemi64alpha")
  test_elecs <- interp_elecs(test_elecs, "A2")
  expect_false(identical(test_elecs$signals$A2, test_data$signals$A2))

})


test_that("Event manipulation works", {

  skip_on_cran()

  expect_is(list_events(test_data), "data.frame")
  test_data <- tag_events(test_data, 255, "osci2")
  test_evs <- list_events(test_data)
  expect_is(test_evs, "data.frame")
  expect_equal(test_data$events$event_label[1], "osci2")
})
