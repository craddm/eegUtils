context("ICA analyses")

demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

demo_SOBI <- run_ICA(demo_epochs)
demo_fastic <- run_ICA(demo_epochs, method = "fastica")

test_that("topoplots for ICA work", {

  skip_on_appveyor()
  vdiffr::expect_doppelganger("topographical plot for SOBI",
                              topoplot(demo_SOBI,
                                       component = "Comp1"))
  vdiffr::expect_doppelganger("topographical plot for fastica",
                              topoplot(demo_fastic,
                                       component = "Comp1"))
})

test_that("ICA timecourses work", {
  vdiffr::expect_doppelganger("timecourse over one component",
                              plot_timecourse(demo_SOBI,
                                              component = "Comp1"))
  vdiffr::expect_doppelganger("fastica timecourse",
                              plot_timecourse(demo_fastic,
                                              component = "Comp2"))
  vdiffr::expect_doppelganger("fastica psd",
                              plot_psd(demo_fastic))
})
