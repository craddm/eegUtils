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
    Fp1 = c(99, 80, 135),
    AFz = c(-299, -200, 35),
    Pz = c(-199, -265, -123),
    time = c( 0, .1, .2)
  )

test_data <- import_raw("Newtest17-256.bdf")

test_that("selection of electrodes and times works as expected", {

  expect_equal(select_elecs(test_dat, electrode = "Fp1"), new_dat)

  expect_equivalent(select_times(test_dat_time, c(-.1, .3)), time_sel)
})

test_that("Selection of electrodes and times works for eeg_* objects", {

  test_Oz <- select_elecs(test_data, electrode = "A2")
  expect_equal(names(test_Oz$signals), "A2")
  expect_s3_class(as.data.frame(test_Oz), "data.frame")

  test_seltime <- select_times(test_data, time_lim = c(-.1, .4))

  expect_true(min(test_seltime$timings$time) >= -.1)
  expect_true(max(test_seltime$timings$time) <= .4)
  expect_equal(nrow(test_seltime$timings), nrow(test_seltime$signals))

  test_epo <- epoch_data(test_data, 255)
  test_epo_elec <- select_elecs(test_epo, electrode = "A5")
  expect_equal(names(test_epo_elec$signals), "A5")
  test_epo_times <- select_times(test_epo, c(-.1, .2))
  min_time <- test_epo$timings$time[which.min(abs(test_epo$timings$time - -.1))]
  max_time <- test_epo$timings$time[which.min(abs(test_epo$timings$time - .2))]
  expect_true(min(test_epo_times$timings$time) >= min_time)
  expect_true(max(test_epo_times$timings$time) <= max_time)

  epo_df <- as.data.frame(test_epo)
  epo_df_A5 <- select_elecs(epo_df, electrode = "A5")
  expect_true("A5" %in% names(epo_df_A5))
  expect_false("A2" %in% names(epo_df_A5))

  test_evoked <- eeg_average(test_epo)
  evo_A5 <- select_elecs(test_evoked, electrode = "A5")
  expect_true("A5" %in% channel_names(evo_A5))
  expect_false("A2" %in% channel_names(evo_A5))

  test_ica <- run_ICA(demo_epochs, pca = 10)
  test_ica_comp1 <- select_elecs(test_ica, "Comp002")
  expect_true("Comp002" %in% channel_names(test_ica_comp1))
  expect_false("Comp001" %in% channel_names(test_ica_comp1))


})

test_that("Selection of epochs functions for eeg_epochs objects only", {

  test_epo <- epoch_data(test_data, 255)
  test_epo <- tag_events(test_epo, 255, "testing")
  test_epo_df <- select_epochs(test_epo, 255, df_out = TRUE)
  expect_is(test_epo_df, "data.frame")
  test_epo_255 <- select_epochs(test_epo, 255)
  expect_is(test_epo_255, "eeg_epochs")
  expect_identical(test_epo_255$signals, test_epo$signals)
  test_epo_testing <- select_epochs(test_epo, "testing")
  expect_identical(test_epo_255$signals,
                   test_epo_testing$signals)
  expect_error(select_epochs(test_epo, 234))
  expect_error(select_epochs(test_epo, "no trig"))
  expect_warning(select_epochs(test_data, 255))
  expect_warning(select_epochs(test_dat, 255))
  sel_epochs <- c(1, 2, 10, 34)
  test_epo_testing <- select_epochs(test_epo,
                                    epoch_no = sel_epochs)
  expect_equal(unique(test_epo_testing$events$epoch),
               sel_epochs)
  expect_equal(unique(test_epo_testing$events$epoch),
               unique(test_epo_testing$timings$epoch))
  test_epo_testing <- select_epochs(test_epo,
                                    epoch_no = sel_epochs,
                                    keep = FALSE)
  expect_equal(unique(test_epo_testing$events$epoch),
               unique(test_epo_testing$timings$epoch))

})
