test_data <- import_raw("Newtest17-256.bdf")
demo_epochs <- electrode_locations(demo_epochs,
                                   montage = "biosemi64alpha",
                                   overwrite = TRUE)
demo_SOBI <- run_ICA(demo_epochs, pca = 10)
demo_tfr <- compute_tfr(demo_epochs,
                        foi = c(4,30),
                        n_freq = 10)

test_that("Plotting of data with multiple epochs works as expected", {
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "epochs plot",
    plot_timecourse(demo_epochs)
    )
  vdiffr::expect_doppelganger(
    "A29 only epochs",
    plot_timecourse(demo_epochs,
                    electrode = "A29")
    )
  vdiffr::expect_doppelganger(
    "A29 baseline corr epochs",
    plot_timecourse(
      demo_epochs,
      baseline = c(-.2, 0),
      electrode = "A29"
    )
  )
  vdiffr::expect_doppelganger(
    "Plot timecourse of component",
    plot_timecourse(demo_SOBI, 2))
  vdiffr::expect_doppelganger(
    "Plot timecourse of evoked",
    plot_timecourse(eeg_average(demo_epochs), 2))
  vdiffr::expect_doppelganger(
    "Plot timecourse of tfr",
    plot_timecourse(eeg_average(demo_tfr)))
  vdiffr::expect_doppelganger(
    "Plot timecourse of tfr in db",
    plot_timecourse(eeg_average(demo_tfr),
                    type = "db",
                    baseline = c(-.1, 0)))
  vdiffr::expect_doppelganger(
    "Plot timecourse of tfr in db at specific freq_range",
    plot_timecourse(eeg_average(demo_tfr),
                    type = "db",
                    baseline = c(-.1, 0),
                    freq_range = c(12, 20)))
})


test_that("Facetting works as expected", {
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "spatial GA plot",
    plot_timecourse(demo_spatial)
  )

  vdiffr::expect_doppelganger(
    "spatial GA plot facetted",
    plot_timecourse(demo_spatial, "Oz", facets = ~epoch_labels)
  )

  vdiffr::expect_doppelganger(
    "spatial GA plot mapping",
    plot_timecourse(demo_spatial, "Oz",
                    mapping = aes(colour = epoch_labels))
  )

  vdiffr::expect_doppelganger(
    "spatial GA plot colour",
    plot_timecourse(demo_spatial, "Oz",
                    colour = "epoch_labels")
  )

  vdiffr::expect_doppelganger(
    "spatial difference plot",
    plot_difference(demo_spatial, electrode = "Oz", conditions = "epoch_labels")
  )
})

test_that("Plotting of butterfly plots from epochs", {
  skip_on_ci()
  vdiffr::expect_doppelganger("butterfly epochs",
                              plot_butterfly(demo_epochs))
  vdiffr::expect_doppelganger("butterfly epochs baseline",
                              plot_butterfly(demo_epochs,
                                             baseline = c(-.2, 0),
                                             electrode = "A29"))
  vdiffr::expect_doppelganger("butterfly evoked",
                              plot_butterfly(eeg_average(demo_epochs)))
})

test_that("erp_raster and erp_image function", {
  skip_on_ci()
  test_epo <- epoch_data(test_data, 255)
  expect_s3_class(erp_image(demo_epochs,
                            electrode = "A29"),
                  "gg")
  expect_s3_class(erp_raster(demo_epochs),
                  "gg")

  expect_error(electrode_locations(test_epo,
                                   montage = "bio2",
                                   overwrite = TRUE))
  test_epo <- electrode_locations(test_epo,
                                  montage = "biosemi64alpha",
                                  overwrite = TRUE)
  expect_s3_class(test_epo$chan_info,
            "data.frame")
})

test_that("erp_scalp runs", {
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "ERP-scalp-plot from demo_spatial",
    erp_scalp(demo_spatial)
  )
  vdiffr::expect_doppelganger(
    "ERP scalp plot with colour",
    erp_scalp(demo_spatial,
              colour = epoch_labels)
  )
})
