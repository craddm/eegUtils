test_that("combining missing part should error", {
  part_001 <- set_participant_id(demo_spatial,
                                 "001")
  part_002 <- set_participant_id(demo_spatial,
                                 "002")
  part_003 <- set_participant_id(demo_spatial,
                                 "003")
  part_003 <- filter(part_003, time < .46)
  part_missing <- demo_spatial
  epochs(part_missing)$participant_id <- NA
  expect_error(eeg_combine(part_001, part_missing))
  expect_warning(eeg_combine(demo_spatial, part_001))
  expect_s3_class(
    eeg_combine(part_002,
                part_001,
                part_003),
    "eeg_group")
  expect_s3_class(eeg_combine(
    list(
      part_002,
      part_001,
      part_003
      )
    ),
    "eeg_group")
  expect_s3_class(eeg_combine(
    list(
      part1 = part_002,
      part2 = part_001,
      part3 = part_003
    )
  ),
  "eeg_group")
})
