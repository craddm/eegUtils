context("Test channel management functions")

temp_chans <- import_chans("standard_1005.elc")

expected_names <- c("electrode",
                    "cart_x",
                    "cart_y",
                    "cart_z",
                    "sph_radius",
                    "sph_phi",
                    "sph_theta",
                    "angle",
                    "radius",
                    "x",
                    "y")
test_that("Import from .elc yields expected file format", {
  expect_equal(expected_names, names(temp_chans))
  expect_equal(cart_to_sph(temp_chans[, 2],
                           temp_chans[, 3],
                           temp_chans[, 4]),
               temp_chans[, 5:7])
  expect_s3_class(plot_electrodes(temp_chans), "ggplot")
  }
)
