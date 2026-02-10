# Plot ERP difference waves

Calculates the difference between the event-related potentials from two
conditions and plots it.

## Usage

``` r
plot_difference(data, ...)

# S3 method for class 'eeg_epochs'
plot_difference(
  data,
  electrode = NULL,
  time_lim = NULL,
  baseline = NULL,
  colour = NULL,
  color = NULL,
  mapping = NULL,
  conditions = "epoch_labels",
  ...
)
```

## Arguments

- data:

  `eegUtils` object. Should have multiple timepoints.

- ...:

  Other arguments passed to methods.

- electrode:

  Electrode(s) to plot.

- time_lim:

  Character vector. Numbers in whatever time unit is used specifying
  beginning and end of time-range to plot. e.g. c(-.1, .3)

- baseline:

  Character vector. Times to use as a baseline. Takes the mean over the
  specified period and subtracts. e.g. c(-.1,0)

- colour:

  Variable to colour lines by. If no variable is passed, only one line
  is drawn.

- color:

  Alias for colour.

- mapping:

  A ggplot2 [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html)
  mapping.

- conditions:

  Defaults to "epoch_labels".

## Value

Returns a `ggplot2` plot object

## Methods (by class)

- `plot_difference(eeg_epochs)`: Plot an ERP difference wave from an
  `eeg_epochs` object

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
plot_difference(demo_spatial, conditions = "epoch_labels", electrode = "P8")
#> Creating epochs based on combinations of variables: participant_id epoch_labels 
#> Ignoring unknown labels:
#> • colour : ""
#> • fill : ""
```
