context("Test import of bdf")

test_data <- import_raw("Newtest17-256.bdf")

test_that("Import and epoching of bdf files works correctly", {

  skip_on_cran()

  expect_s3_class(test_data, "eeg_data")
  expect_true(is.eeg_data(test_data))
  expect_false(is.eeg_epochs(test_data))
  expect_true(is.data.frame(as.data.frame(test_data)))

  test_epo <- epoch_data(test_data, 255)
  expect_s3_class(test_epo, "eeg_epochs")

  test_evo <- eeg_average(test_epo)
  expect_s3_class(test_evo, "eeg_evoked")

})
