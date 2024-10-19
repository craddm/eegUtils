# Helper function to create dummy eeg_data objects
create_dummy_eeg_data <- function(participant_id,
                                  n_samples = 100,
                                  n_channels = 3) {
  data <- list(
    signals = data.frame(matrix(rnorm(
      n_samples * n_channels
    ), ncol = n_channels)),
    events = data.frame(
      event_onset = 1:5,
      event_time = (1:5) / 100,
      event_type = 11:15
    ),
    timings = data.frame(
      sample = 1:n_samples,
      time = (1:n_samples - 1) / 100
    ),
    chan_info = data.frame(channel = paste0("chan", 1:n_channels)),
    srate = 100
  )
  class(data) <- c("eeg_data")
  data$epochs <- data.frame(participant_id = participant_id,
                            epoch = 1,
                            recording = "test")
  data
}

# Helper function to create dummy eeg_epochs objects
create_dummy_eeg_epochs <- function(participant_id,
                                    n_epochs = 5,
                                    n_samples = 100,
                                    n_channels = 3) {
  data <- list(
    signals = data.frame(matrix(
      rnorm(n_epochs * n_samples * n_channels), ncol = n_channels
    )),
    events = data.frame(
      event_onset = rep(1:n_samples, n_epochs),
      event_time = rep((1:n_samples - 1) / 100, n_epochs),
      epoch = rep(1:n_epochs, each = n_samples),
      time = rep(0, n_epochs)
    ),
    timings = data.frame(
      sample = 1:n_samples,
      time = (1:n_samples - 1) / 100,
      epoch = rep(1:n_epochs, each = n_samples)
    ),
    chan_info = data.frame(channel = paste0("chan", 1:n_channels)),
    epochs = data.frame(
      participant_id = participant_id,
      epoch = 1:n_epochs,
      recording = "test"
    ),
    srate = 100
  )
  class(data) <- c("eeg_epochs", "eeg_data")
  data
}

test_that("eeg_combine works with eeg_data objects", {
  data1 <- create_dummy_eeg_data("P1")
  data2 <- create_dummy_eeg_data("P1")

  combined <- eeg_combine(data1, data2)

  expect_s3_class(combined, "eeg_data")
  expect_equal(nrow(combined$signals),
               nrow(data1$signals) + nrow(data2$signals))
  expect_equal(nrow(combined$events),
               nrow(data1$events) + nrow(data2$events))
  expect_equal(nrow(combined$timings),
               nrow(data1$timings) + nrow(data2$timings))
})

test_that("eeg_combine works with eeg_epochs objects", {
  data1 <- create_dummy_eeg_epochs("P1")
  data2 <- create_dummy_eeg_epochs("P1")

  combined <- eeg_combine(data1, data2)

  expect_s3_class(combined, "eeg_epochs")
  expect_equal(nrow(combined$signals),
               nrow(data1$signals) + nrow(data2$signals))
  expect_equal(nrow(combined$events),
               nrow(data1$events) + nrow(data2$events))
  expect_equal(nrow(combined$epochs),
               nrow(data1$epochs) + nrow(data2$epochs))
})

test_that("eeg_combine creates eeg_group for multiple participants", {
  data1 <- create_dummy_eeg_epochs("P1")
  data2 <- create_dummy_eeg_epochs("P2")

  combined <- eeg_combine(data1, data2)

  expect_s3_class(combined, "eeg_group")
  expect_equal(length(unique(combined$epochs$participant_id)), 2)
})

test_that("eeg_combine checks for participant_id", {
  data1 <- create_dummy_eeg_epochs("P1")
  data2 <- create_dummy_eeg_epochs("P2")
  data2$epochs$participant_id <- NA

  expect_error(eeg_combine(data1, data2), "`participant_id` is missing")
})

test_that("check_timings corrects epoch numbering", {
  data1 <- create_dummy_eeg_epochs("P1", n_epochs = 3)
  data2 <- create_dummy_eeg_epochs("P1", n_epochs = 3)
  data2$epochs$epoch <- 1:3
  data2$timings$epoch <- rep(1:3, each = 100)

  combined <- eeg_combine(data1, data2, check_timings = TRUE)

  expect_equal(unique(combined$epochs$epoch), 1:6)
  expect_equal(max(combined$timings$epoch), 6)
})

test_that("eeg_combine.list works with a list of eeg objects", {
  data_list <- list(
    create_dummy_eeg_epochs("P1"),
    create_dummy_eeg_epochs("P2"),
    create_dummy_eeg_epochs("P3")
  )

  combined <- eeg_combine(data_list)

  expect_s3_class(combined, "eeg_group")
  expect_equal(length(unique(combined$epochs$participant_id)), 3)
})

test_that("eeg_combine throws error for incompatible objects", {
  data1 <- create_dummy_eeg_data("P1")
  data2 <- create_dummy_eeg_epochs("P2")

  expect_error(eeg_combine(data1, data2),
               "All objects must be unepoched eeg_data objects")
})
