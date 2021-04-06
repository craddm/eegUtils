test_that("import_set files", {
  check_new_eeglab <- import_set("small_test_hdf5.set")
  expect_s3_class(
    check_new_eeglab,
    "eeg_epochs"
    )
})
