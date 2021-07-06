test_that("epoch images plot correctly", {
  skip_on_ci()
  vdiffr::expect_doppelganger(
    "ERP image A29 from demo epochs",
    erp_image(demo_epochs,
              "A29")
  )
})
