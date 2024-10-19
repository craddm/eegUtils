test_that("glm fitting works", {
  test_glm <- readRDS("reference_files/test_glm.rds")
  expect_equal(fit_glm(~epoch_labels,
                       data = demo_spatial),
               readRDS("reference_files/test_glm.rds"),
               ignore_attr = TRUE)
  expect_equal(fit_glm(~epoch_labels + baseline,
                       data = demo_spatial,
                       baseline = c(-.1, 0)),
               readRDS("reference_files/test_glm_baseline.rds"),
               ignore_attr = TRUE)
})
