context ("Test selection")

test_dat <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555),
    AFz = c(120,-299,-200, 35, 2112,-3, 243),
    Pz = c(-24,-199,-265,-123, 12, 103, 143)
  )

test_dat_time <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555),
    AFz = c(120,-299,-200, 35, 2112,-3, 243),
    Pz = c(-24,-199,-265,-123, 12, 103, 143),
    time = c(-.1, 0, .1, .2, .3, .4, .5)
  )

new_dat <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555)
  )

time_sel <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515),
    AFz = c(120,-299,-200, 35, 2112),
    Pz = c(-24,-199,-265,-123, 12),
    time = c(-.1, 0, .1, .2, .3)
  )

test_that("selection of electrodes and times works as expected", {

  expect_equal(select_elecs(test_dat, electrode = "Fp1"), new_dat)

  expect_equal(select_times(test_dat_time, c(-.1, .3)), time_sel)
})
