test_that("epoch images plot correctly", {
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "ERP image A29 from demo epochs",
    erp_image(demo_epochs,
              "A29"))
  vdiffr::expect_doppelganger(
    "ERP raster with facets",
    erp_raster(demo_spatial) + facet_wrap(~epoch_labels)
  )
  vdiffr::expect_doppelganger(
    "ERP raster without anat_order",
    erp_raster(demo_spatial, anat_order = FALSE) + facet_wrap(~epoch_labels)
  )
})
