context("Test artefact rejection methods")

test_that("FASTER runs correctly.", {
  test_FASTER <- ar_FASTER(demo_epochs)
  expect_known_output(test_FASTER, "faster_t.Rdata")
  #expect_known_hash(test_FASTER, "04de56b5ba")
})

test_that("Calculate epoch stats", {
  demo_epo <- epoch_stats(demo_epochs)
  expect_known_output(demo_epo, "epo_stats.Rdata")
})

test_that("Calculating channel stats", {
  chan_stats <- channel_stats(demo_epochs)
  expect_known_output(chan_stats, "chan_stats.Rdata")
})

test_that("ar_thresh runs correctly.", {
  test_thresh <- ar_thresh(demo_epochs,
                           50)
  expect_known_output(test_thresh,
                      "thresh_test.Rdata")
})
