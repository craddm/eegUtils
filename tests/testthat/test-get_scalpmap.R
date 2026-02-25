test_that("point_in_hull classifies points correctly", {
  # square hull: (1,1), (1,-1), (-1,-1), (-1,1)
  elec_x <- c(1, 1, -1, -1)
  elec_y <- c(1, -1, -1, 1)

  # inside point
  expect_true(point_in_hull(elec_x, elec_y, 0, 0))
  # outside point
  expect_false(point_in_hull(elec_x, elec_y, 2, 0))
  # another outside point

  expect_false(point_in_hull(elec_x, elec_y, 0, 2))

  # vectorized: mix of inside and outside
  gx <- c(0, 0.5, 2, -0.5, 0, 1.5)
  gy <- c(0, 0.5, 0, -0.5, 1.5, 1.5)
  result <- point_in_hull(elec_x, elec_y, gx, gy)
  expect_equal(result, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
})

test_that("biharmonic respects convex hull boundary", {
  # Create synthetic data with known square electrode positions
  synth <- data.frame(
    x = c(1, 1, -1, -1, 0),
    y = c(1, -1, -1, 1, 0),
    fill = c(1, 2, 3, 4, 2.5),
    electrode = c("A", "B", "C", "D", "E")
  )

  result <- biharmonic(synth,
                       grid_res = 20,
                       interp_limit = "convex_hull",
                       r = 3)

  # All returned points should be inside the convex hull (with buffer)
  in_hull <- point_in_hull(synth$x, synth$y, result$x, result$y)
  expect_true(all(in_hull))

  # A point like (1.5, 0) should not appear in the result
  outside_pts <- result[result$x > 1.1 & abs(result$y) < 0.1, ]
  expect_equal(nrow(outside_pts), 0)
})

test_that("convex_hull produces subset of head region", {
  skip_on_cran()
  demo_epochs <- electrode_locations(demo_epochs,
                                     montage = "biosemi64alpha",
                                     overwrite = TRUE)

  head_result <- get_scalpmap(demo_epochs,
                              interp_limit = "head",
                              grid_res = 50)
  hull_result <- get_scalpmap(demo_epochs,
                              interp_limit = "convex_hull",
                              grid_res = 50)

  expect_lte(nrow(hull_result), nrow(head_result))
})

test_that("fit_gam_topo works with convex_hull", {
  skip_on_cran()
  demo_epochs <- electrode_locations(demo_epochs,
                                     montage = "biosemi64alpha",
                                     overwrite = TRUE)

  demo_df <- as.data.frame(demo_epochs, long = TRUE)
  demo_df <- dplyr::summarise(
    dplyr::group_by(demo_df, x, y, electrode),
    fill = mean(amplitude, na.rm = TRUE),
    .groups = "drop"
  )
  demo_df <- demo_df[stats::complete.cases(demo_df), ]

  result <- fit_gam_topo(demo_df,
                         grid_res = 50,
                         interp_limit = "convex_hull",
                         r = 95,
                         k = 10)

  in_hull <- point_in_hull(demo_df$x, demo_df$y, result$x, result$y)
  expect_true(all(in_hull))
})

test_that("get_scalpmap accepts convex_hull", {
  skip_on_cran()
  demo_epochs <- electrode_locations(demo_epochs,
                                     montage = "biosemi64alpha",
                                     overwrite = TRUE)

  result <- get_scalpmap(demo_epochs,
                         interp_limit = "convex_hull",
                         grid_res = 50)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(c("x", "y", "fill") %in% names(result)))
})
