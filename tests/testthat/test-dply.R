context("Test dplyr exports")

test_that("selection of electrodes and times works as expected", {

  add_col <- function(.data) {
    .data$signals$yoyo <-
      (.data$signals$A5 + .data$signals$A13 + .data$signals$A29) / 3
    .data
    }

  expect_equivalent(select(demo_epochs, A5),
                    select_elecs(demo_epochs, "A5"))

  demo_epochs$signals <- tibble::as_tibble(demo_epochs$signals)
  rownames(demo_epochs$signals) <- NULL
  expect_equivalent(filter(demo_epochs, time >= -.1, time <= .3),
               select_times(demo_epochs, c(-.1, .3)))
  expect_equal(mutate(demo_epochs, yoyo = (A5 + A13 + A29) / 3),
               add_col(demo_epochs))
  expect_equivalent(filter(demo_epochs, epoch <= 10, epoch >= 5),
                    select_epochs(demo_epochs, epoch_no = 5:10))
})
