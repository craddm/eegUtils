context("Test channel management functions")

temp_chans <- import_chans("standard_1005.elc")

expected_names <- c("electrode",
                    "radius",
                    "theta",
                    "phi",
                    "cart_x",
                    "cart_y",
                    "cart_z",
                    "x",
                    "y")
test_that("Import from .elc yields expected file format", {
  expect_equal(expected_names, names(temp_chans))
  cts <- cart_to_spherical(temp_chans[, c("cart_x", "cart_y", "cart_z")])
  #cts[, 2:3] <- cts[, 2:3] * 180 / pi
  # expect_equal(tibble::as.tibble(cts),
  #              temp_chans[, c("radius",
  #                             "theta",
  #                             "phi")])
  expect_s3_class(plot_electrodes(temp_chans), "ggplot")
  }
)
