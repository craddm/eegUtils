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

test_that("ar_peak2peak marks epochs correctly", {
  # threshold = 30 catches epochs since demo max p2p ~ 64 uV
  result <- ar_peak2peak(demo_epochs, threshold = 30)
  expect_s3_class(result, "eeg_epochs")
  expect_true(!is.null(result$reject$epochs))
  expect_true(all(result$reject$epochs$reason == "peak-to-peak"))
})

test_that("ar_peak2peak rejects epochs when reject = TRUE", {
  result <- ar_peak2peak(demo_epochs, threshold = 30, reject = TRUE)
  expect_s3_class(result, "eeg_epochs")
  original_epochs <- length(unique(demo_epochs$timings$epoch))
  result_epochs <- length(unique(result$timings$epoch))
  expect_lt(result_epochs, original_epochs)
})

test_that("ar_peak2peak works with a high threshold (no epochs rejected)", {
  result <- ar_peak2peak(demo_epochs, threshold = 1e6)
  expect_s3_class(result, "eeg_epochs")
  expect_null(result$reject$epochs)
})

test_that("ar_peak2peak snapshot", {
  expect_snapshot(ar_peak2peak(demo_epochs, threshold = 30))
})

test_that("ar_gradient marks epochs correctly", {
  # threshold = 10 catches epochs since demo max gradient ~ 19.9 uV
  result <- ar_gradient(demo_epochs, threshold = 10)
  expect_s3_class(result, "eeg_epochs")
  expect_true(!is.null(result$reject$epochs))
  expect_true(all(result$reject$epochs$reason == "gradient"))
})

test_that("ar_gradient rejects epochs when reject = TRUE", {
  result <- ar_gradient(demo_epochs, threshold = 10, reject = TRUE)
  expect_s3_class(result, "eeg_epochs")
  original_epochs <- length(unique(demo_epochs$timings$epoch))
  result_epochs <- length(unique(result$timings$epoch))
  expect_lt(result_epochs, original_epochs)
})

test_that("ar_gradient works with a high threshold (no epochs rejected)", {
  result <- ar_gradient(demo_epochs, threshold = 1e6)
  expect_s3_class(result, "eeg_epochs")
  expect_null(result$reject$epochs)
})

test_that("ar_gradient snapshot", {
  expect_snapshot(ar_gradient(demo_epochs, threshold = 10))
})

test_that("ar_flat marks epochs correctly", {
  # threshold = 10 catches epochs since demo variance min ~ 4.1 uV^2
  result <- ar_flat(demo_epochs, threshold = 10)
  expect_s3_class(result, "eeg_epochs")
  expect_true(!is.null(result$reject$epochs))
  expect_true(all(result$reject$epochs$reason == "flat"))
})

test_that("ar_flat rejects epochs when reject = TRUE", {
  result <- ar_flat(demo_epochs, threshold = 10, reject = TRUE)
  expect_s3_class(result, "eeg_epochs")
  original_epochs <- length(unique(demo_epochs$timings$epoch))
  result_epochs <- length(unique(result$timings$epoch))
  expect_lt(result_epochs, original_epochs)
})

test_that("ar_flat works with a very low threshold (no epochs rejected)", {
  result <- ar_flat(demo_epochs, threshold = 0)
  expect_s3_class(result, "eeg_epochs")
  expect_null(result$reject$epochs)
})

test_that("ar_flat snapshot", {
  expect_snapshot(ar_flat(demo_epochs, threshold = 10))
})

test_that("ar_* functions error on unsupported class", {
  expect_error(ar_peak2peak(list(), threshold = 100))
  expect_error(ar_gradient(list(), threshold = 50))
  expect_error(ar_flat(list(), threshold = 0.01))
})
