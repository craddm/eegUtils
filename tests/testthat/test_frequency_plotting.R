test_that("plot_psd produces a ggplot", {
  psd_test <- plot_psd(demo_epochs)
  vdiffr::expect_doppelganger("power spectrum for demo epochs",
                              psd_test)
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 16,
                          n_cycles = 3,
                          output = "fourier")
  tfr_plot <- plot_tfr(tfr_test)
  vdiffr::expect_doppelganger("TFA plot for demo epochs fourier",
                              tfr_plot)
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 16,
                          n_cycles = 3)
  tfr_plot <- plot_tfr(tfr_test)
  vdiffr::expect_doppelganger("TFA plot for demo epochs",
                              tfr_plot)
  baseline_types <- c("db", "ratio", "pc", "absolute")
  tfr_plots <- purrr::map(baseline_types,
                          ~plot_tfr(tfr_test,
                                    baseline_type = .,
                                    baseline = c(-.1, 0)))
  purrr::walk(1:length(baseline_types),
             ~vdiffr::expect_doppelganger(baseline_types[.],
                                          tfr_plots[.]))
  #vdiffr::expect_doppelganger("TFA plot for demo epochs", tfr_plot)
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 16,
                          n_cycles = 3,
                          keep_trials = TRUE)
  tfr_plot <- plot_tfr(tfr_test)
  vdiffr::expect_doppelganger("TFA plot unaveraged", tfr_plot)
  tfr_test <- compute_tfr(demo_epochs,
                          foi = c(4, 30),
                          n_freq = 16,
                          n_cycles = 3, spacing = "log")
  tfr_plot <- plot_tfr(tfr_test)
  vdiffr::expect_doppelganger("TFA plot for log spaced",
                              tfr_plot)
})

tfr_test <- compute_tfr(demo_epochs,
                        foi = c(4, 30),
                        n_freq = 16,
                        n_cycles = 3)

test_that("Timecourses from eeg_tfr objects work", {
  vdiffr::expect_doppelganger("timecourses from `eeg_tfr`",
                              plot_timecourse(tfr_test))
  vdiffr::expect_doppelganger("timecourses from `eeg_tfr` with freq_range",
                              plot_timecourse(tfr_test,
                                              freq_range = c(8, 12)))
  vdiffr::expect_doppelganger("timecourses from `eeg_tfr` with time_lim",
                              plot_timecourse(tfr_test,
                                              freq_range = c(8, 12),
                                              time_lim = c(-.2, .5)))
  vdiffr::expect_doppelganger("timecourses from `eeg_tfr` with baseline",
                              plot_timecourse(tfr_test,
                                              freq_range = c(15, 30),
                                              baseline = c(-.3, -.1),
                                              time_lim = c(-.2, .5)))
})
