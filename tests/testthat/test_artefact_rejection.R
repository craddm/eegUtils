test_data <- import_raw("Newtest17-256.bdf")

test_that("FASTER runs correctly.", {
  test_FASTER <- ar_FASTER(demo_epochs)
  expect_equal(test_FASTER,
               readRDS("reference_files/test_faster.rds"),
               ignore_attr = TRUE)
})

test_that("Calculate epoch stats", {
  demo_epo <- epoch_stats(demo_epochs)
  expect_snapshot(demo_epo)
})

test_that("Calculating channel stats", {
  chan_stats <- channel_stats(demo_epochs)
  expect_snapshot(chan_stats)
})

test_that("ar_thresh runs correctly.", {
  expect_snapshot(ar_thresh(demo_epochs,
                            20))
  expect_equal(ar_thresh(demo_epochs,
                         20),
               readRDS("reference_files/test_thresh.rds"),
               ignore_attr = TRUE)
})

test_that("ar_thresh runs correctly.", {
  expect_snapshot(ar_thresh(test_data,
                            30))
})
