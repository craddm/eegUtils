# Time-frequency plot

Creates a time-frequency plot of an `eeg_tfr` object. The plot has time
on the x-axis and frequency on the y-axis. If no electrode is supplied,
it will average over all electrodes.

## Usage

``` r
plot_tfr(
  data,
  electrode = NULL,
  time_lim = NULL,
  freq_range = NULL,
  baseline_type = "none",
  baseline = NULL,
  fill_lims = NULL,
  interpolate = FALSE,
  na.rm = TRUE,
  fun.data = mean
)
```

## Arguments

- data:

  Object of class `eeg_tfr`

- electrode:

  Electrode to plot. If none is supplied, averages over all electrodes.

- time_lim:

  Time limits of plot.

- freq_range:

  Vector of two numbers. (e.g. c(8, 40)).

- baseline_type:

  baseline correction to apply. Defaults to `none`.

- baseline:

  Baseline period

- fill_lims:

  Custom colour scale (i.e. range of power). e.g. c(-5, 5).

- interpolate:

  Interpolation of raster for smoother plotting.

- na.rm:

  Remove NA values silently (TRUE) or with a warning (FALSE). Defaults
  to TRUE.

- fun.data:

  Statistical function to use for averaging over electrodes/conditions.
  Defaults to `mean`. Alternatively, supply `weighted.mean`, which will
  attempt to weight separate conditions appropriately.

## Value

A `ggplot` object

## Details

Various different baseline options can be applied here (e.g. "db" for
decibels, "pc" for percent change, "divide" for division; see
`rm_baseline` for details).

## See also

[`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)

## Author

Matt Craddock <matt@mattcraddock.com>
