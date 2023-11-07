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
  demo_reref <- eeg_reference(demo_epochs, "A29")
  expect_false("A29" %in% names(demo_reref$signals))
  expect_true(identical("A29", demo_reref$reference$ref_chans))
})
