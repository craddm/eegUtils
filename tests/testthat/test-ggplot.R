context("Plotting functions")

demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

test_that("geom_topo testing", {
  skip_on_ci()
  vdiffr::expect_doppelganger("geom_topo_test",
                              ggplot(demo_epochs,
                                     aes(x = x,
                                         y = y,
                                         fill = amplitude,
                                         z = amplitude)) +
                                geom_topo())
  vdiffr::expect_doppelganger("geom_topo_name",
                              ggplot(demo_epochs,
                                     aes(x = x,
                                         y = y,
                                         fill = amplitude,
                                         z = amplitude)) +
                                geom_topo(chan_markers = "text",
                                          aes(label = electrode)))
  vdiffr::expect_doppelganger("geom_topo_build",
                              ggplot(demo_epochs,
                                     aes(x = x,
                                         y = y,
                                         fill = amplitude)) +
                                stat_scalpmap() +
                                geom_mask(r = 115) +
                                geom_head() +
                                geom_channels())
  vdiffr::expect_doppelganger("geom_topo_labels",
                              ggplot(demo_epochs,
                                     aes(x = x,
                                         y = y,
                                         fill = amplitude)) +
                                stat_scalpmap() +
                                geom_mask(r = 115) +
                                geom_head() +
                                geom_channels(geom = "text",
                                              aes(label = electrode),
                                              size = 6)
                              )
  vdiffr::expect_doppelganger("geom_topo_head_test",
                              ggplot(demo_epochs,
                                     aes(x = x,
                                         y = y,
                                         fill = amplitude,
                                         z = amplitude)) +
                                geom_topo(interp_limit = "head"))
  vdiffr::expect_doppelganger("geom_statscalp_head",
                              ggplot(demo_epochs,
                                     aes(x = x,
                                         y = y,
                                         fill = amplitude)) +
                                stat_scalpmap(interp_limit = "head") +
                                geom_mask(r = 94) +
                                geom_head(r = 92) +
                                geom_channels(geom = "text",
                                              aes(label = electrode))
  )
  }
)
