context("Test channel management functions")

temp_chans <- import_chans("standard_1005.elc")

expected_names <- c("electrode",
                    "cart_x",
                    "cart_y",
                    "cart_z",
                    "sph_radius",
                    "sph_phi",
                    "sph_theta",
                    "pol_theta",
                    "pol_radius",
                    "angle",
                    "radius",
                    "x",
                    "y")
test_that("Import from .elc yields expected file format", {
  expect_equal(expected_names, names(temp_chans))
  cts <- cart_to_sph(temp_chans[, 2, drop = TRUE],
              temp_chans[, 3, drop = TRUE],
              temp_chans[, 4, drop = TRUE])
  #cts[, 2:3] <- cts[, 2:3] * 180 / pi
  expect_equal(tibble::as.tibble(cts),
               temp_chans[, 5:7])
  expect_s3_class(plot_electrodes(temp_chans), "ggplot")
  }
)
