context("ICA analyses")

demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

test_that("SOBI ICA can be calculated and plotted", {
  demo_SOBI <- run_ICA(demo_epochs)
  vdiffr::expect_doppelganger("topographical plot for SOBI",
                              topoplot(demo_SOBI,
                                       component = "Comp1"))
  vdiffr::expect_doppelganger("timecourse over one component",
                              plot_timecourse(demo_SOBI,
                                              component = "Comp1"))
})
