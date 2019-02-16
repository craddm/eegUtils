context("Test data modification")

test_data <- import_raw("Newtest17-256.bdf")

test_that("Referencing works for eeg_* objects", {

  skip_on_cran()

  ref_test <- rowMeans(test_data$signals)
  test_reref <- eeg_reference(test_data)

  expect_equal(as.data.frame(test_data$signals),
               test_reref$signals + ref_test)
  expect_true(all(c("ref_chans", "excluded") %in% names(test_reref$reference)))

  test_reref <- eeg_reference(test_data, exclude = "A15")
  expect_match(test_reref$reference$excluded, "A15")
  test_reref <- eeg_reference(test_data, ref_chans = "A2")
  expect_false("A2" %in% names(test_reref$signals))
  test_reref <- eeg_reference(test_data, "A14")
  expect_true("A2" %in% names(test_reref$signals))
  expect_false("A14" %in% names(test_reref$signals))
})


test_that("Removing baseline works", {

  skip_on_cran()

  test_bl <- rm_baseline(test_data)
  expect_equal(test_bl$signals$A1,
               test_data$signals$A1 - mean(test_data$signals$A1))
  test_epo <- epoch_data(test_data, 255)
  test_bl <- rm_baseline(test_epo, c(-.1, 0))

})

test_that("Downsampling output is sensible", {

  test_epo <- epoch_data(test_data, 255)
  test_reref <- eeg_reference(test_epo)
  test_ds <- eeg_downsample(test_reref, 2)
  expect_equal(nrow(test_ds$signals), nrow(test_ds$timings))
})
