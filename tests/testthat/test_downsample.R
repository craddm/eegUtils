context("Test downsampling")

test_data <- import_raw("Newtest17-256.bdf")

test_that("Downsampling output is sensible", {

  skip_on_cran()

  test_epo <- epoch_data(test_data, 255)
  test_ds <- eeg_downsample(test_epo, 2)
  expect_equal(nrow(test_ds$signals), nrow(test_ds$timings))
})


