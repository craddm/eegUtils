context ("Test selection")

test_dat <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555),
    AFz = c(120, -299, -200, 35, 2112, -3, 243),
    Pz = c(-24, -199, -265, -123, 12, 103, 143)
  )

test_dat_time <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555),
    AFz = c(120, -299, -200, 35, 2112, -3, 243),
    Pz = c(-24, -199, -265, -123, 12, 103, 143),
    time = c(-.1, 0, .1, .2, .3, .4, .5)
  )

new_dat <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515, 3, 555)
  )

time_sel <-
  data.frame(
    Fp1 = c(100, 99, 80, 135, 3515),
    AFz = c(120, -299, -200, 35, 2112),
    Pz = c(-24, -199, -265, -123, 12),
    time = c(-.1, 0, .1, .2, .3)
  )
test_data <- import_raw("Newtest17-256.bdf")

test_that("selection of electrodes and times works as expected", {

  expect_equal(select_elecs(test_dat, electrode = "Fp1"), new_dat)

  expect_equal(select_times(test_dat_time, c(-.1, .3)), time_sel)
})

test_that("Selection of electrodes and times works for eeg_* objects",{

  test_Oz <- select_elecs(test_data, electrode = "A2")
  expect_equal(names(test_Oz$signals), "A2")
  expect_s3_class(as.data.frame(test_Oz), "data.frame")

  test_seltime <- select_times(test_data, time_lim = c(-.1, .2))

  expect_true(min(test_seltime$timings$time) >= -.1)
  expect_true(max(test_seltime$timings$time) <= .2)
  expect_equal(nrow(test_seltime$timings), nrow(test_seltime$signals))

  test_epo <- epoch_data(test_data, 255)
  test_epo_elec <- select_elecs(test_epo, electrode = "A5")
  expect_equal(names(test_epo_elec$signals), "A5")
  test_epo_times <- select_times(test_epo, c(-.1, .2))
  min_time <- test_epo$timings$time[which.min(abs(test_epo$timings$time - -.1))]
  max_time <- test_epo$timings$time[which.min(abs(test_epo$timings$time - .2))]
  expect_true(min(test_epo_times$timings$time) >= min_time)
  expect_true(max(test_epo_times$timings$time) <= max_time)

})
