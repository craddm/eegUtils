test_data <- import_raw("Newtest17-256.bdf")
test_that("Removing baseline works", {

  skip_on_cran()

  test_bl <- rm_baseline(test_data)
  expect_equal(test_bl$signals$A1,
               test_data$signals$A1 - mean(test_data$signals$A1))
  test_epo <- epoch_data(test_data, 255)
  test_bl <- rm_baseline(test_epo,
                         c(-.1, 0))
  expect_equal(test_bl, readRDS("reference_files/baselined_test.rds"))
})
