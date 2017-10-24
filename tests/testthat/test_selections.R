context ("Test selection")

test_dat <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555),
    AFz = c(120,-299,-200, 35, 2112,-3, 243),
    Pz = c(-24,-199,-265,-123, 12, 103, 143)
  )

new_dat <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555)
  )

test_that("selection of electrodes works as expected", {
  expect_equal(select_elecs(test_dat, electrode = "Fp1"), new_dat)
})
