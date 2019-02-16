context("Plotting functions")

demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

test_that("geom_topo testing", {
          vdiffr::expect_doppelganger("geom_topo_test",
                                      ggplot(demo_epochs,
                                             aes(x = x,
                                                 y = y,
                                                 fill = amplitude)) +
                                        geom_topo())
          vdiffr::expect_doppelganger("geom_topo_build",
                                      ggplot(demo_epochs,
                                             aes(x = x,
                                                 y = y,
                                                 fill = amplitude)) +
                                        stat_scalpmap() +
                                        geom_mask(scale_fac = 1.4) +
                                        geom_head() +
                                        geom_channels())
          }
)
