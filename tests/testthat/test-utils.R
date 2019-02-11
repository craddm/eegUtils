context("Test utils.R")

cond_labs_l <- "Light"
cond_labs_nl <- "No light"
cond_labs_lsf <- "LSF"
cond_labs_l_lsf <- "Light/LSF"
cond_labs_sf <- c("HSF", "LSF")
cond_labs_mix <- c("HSF", "Match/LSF")

orig_labs <- c("Light/HSF", "Light/LSF", "No light/HSF", "No light/LSF")


test_that("Label checking", {
  expect_true(label_check(cond_labs_l, orig_labs))
  expect_true(label_check(cond_labs_nl, orig_labs))
  expect_true(label_check(cond_labs_lsf, orig_labs))
  expect_true(label_check(cond_labs_l_lsf, orig_labs))
  expect_true(all(label_check(cond_labs_sf, orig_labs)))
  expect_error(label_check(cond_labs_mix, orig_labs))
})

test_that("Padding and unpadding are functional", {
  x <- rnorm(10)
  expect_equal(length(pad(x, 10)), length(x) + 20)
  expect_equal(unpad(pad(x, 10), 10), x)
})

test_that("Sample n calculated correctly", {
  expect_equal(samples(demo_epochs),
               nrow(demo_epochs$signals))
})
