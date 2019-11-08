context("Test filtering")

test_data <- import_raw("Newtest17-256.bdf")

test_that("filtering behaves", {

  skip_on_cran()

  test_epo <- epoch_data(test_data, 255)
  test_ds <- eeg_downsample(test_epo, 2)
  expect_equal(nrow(test_ds$signals), nrow(test_ds$timings))
})


