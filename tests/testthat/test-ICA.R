context("ICA analyses")

demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)

demo_SOBI <- run_ICA(demo_epochs)
demo_fastic <- run_ICA(demo_epochs, method = "fastica")
demo_fica <- run_ICA(demo_epochs, method = "fica", pca = 10)

test_that("topoplots for ICA work", {

  skip_on_appveyor()
  vdiffr::expect_doppelganger("topographical plot for SOBI",
                              topoplot(demo_SOBI,
                                       component = "Comp1"))
  vdiffr::expect_doppelganger("topographical plot for fastica",
                              topoplot(demo_fastic,
                                       component = "Comp1"))
  vdiffr::expect_doppelganger("topographical plot for fica",
                              topoplot(demo_fica,
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

test_that("component removal works", {
  demo_rm <- apply_ica(demo_epochs, demo_SOBI, 1)
  demo_re <- apply_ica(demo_SOBI)
  vdiffr::expect_doppelganger("removed one component",
                              plot_butterfly(demo_rm))
  vdiffr::expect_doppelganger("reconstruct all",
                              plot_butterfly(demo_re))
})
