demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)
set.seed(1000)
demo_SOBI <- run_ICA(demo_epochs, pca = 10)
demo_fastic <- run_ICA(demo_epochs, method = "fastica", pca = 10)
demo_fica <- run_ICA(demo_spatial, method = "fica", pca = 14)

test_that("topoplots for ICA work", {

  skip_on_ci()
  vdiffr::expect_doppelganger("topographical plot for SOBI",
                              topoplot(demo_SOBI,
                                       component = "Comp001",
                                       limits = c(-4, 5), verbose = FALSE))
  vdiffr::expect_doppelganger("topographical plot for fastica",
                              topoplot(demo_fastic,
                                       component = "Comp001",
                                       limits = c(-5, 5)))
  vdiffr::expect_doppelganger("topographical plot for fica",
                              topoplot(demo_fica,
                                       component = "Comp001",
                                       limits = c(-4, 3.5)))
  vdiffr::expect_doppelganger("topoplot multiple comps",
                              topoplot(demo_fica,
                                       component = c("Comp001", "Comp002"),
                                       limits = c(-4, 3.5)))
})

test_that("ICA timecourses work", {

  vdiffr::expect_doppelganger("timecourse over one component",
                              plot_timecourse(demo_SOBI,
                                              component = "Comp001"))
  vdiffr::expect_doppelganger("fastica timecourse",
                              plot_timecourse(demo_fastic,
                                              component = "Comp002"))
  vdiffr::expect_doppelganger("fastica psd",
                              plot_psd(demo_fastic, verbose = FALSE))
})

test_that("component removal works", {
  demo_rm <- apply_ica(demo_epochs, demo_SOBI, 1)
  demo_re <- apply_ica(demo_SOBI)
  vdiffr::expect_doppelganger("removed one component",
                              plot_butterfly(demo_rm))
  vdiffr::expect_doppelganger("reconstruct all",
                              plot_butterfly(demo_re))
})

test_that("artefact detect works", {

  expect_equal(ar_acf(demo_SOBI, plot = FALSE), character())
  expect_equal(ar_chanfoc(demo_SOBI, plot = FALSE), character())
  expect_equal(ar_trialfoc(demo_SOBI, plot = FALSE), "Comp009")
  expect_equal(ar_eogcor(demo_fica,
                         demo_spatial,
                         HEOG = c("EXG1", "EXG2"),
                         VEOG = c("EXG3", "EXG4"),
                         plot = FALSE),
               character(0))

})
