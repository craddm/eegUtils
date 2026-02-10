# Create a butterfly plot from timecourse data

Typically event-related potentials/fields, but could also be timecourses
from frequency analyses for single frequencies. Output is a ggplot2
object. CIs not possible.

## Usage

``` r
plot_butterfly(data, ...)

# Default S3 method
plot_butterfly(
  data,
  time_lim = NULL,
  baseline = NULL,
  legend = TRUE,
  continuous = FALSE,
  browse_mode = FALSE,
  allow_facets = FALSE,
  ...
)

# S3 method for class 'eeg_evoked'
plot_butterfly(
  data,
  time_lim = NULL,
  baseline = NULL,
  legend = TRUE,
  continuous = FALSE,
  browse_mode = FALSE,
  allow_facets = FALSE,
  ...
)

# S3 method for class 'eeg_data'
plot_butterfly(
  data,
  time_lim = NULL,
  baseline = NULL,
  legend = TRUE,
  allow_facets = FALSE,
  browse_mode = FALSE,
  ...
)

# S3 method for class 'eeg_epochs'
plot_butterfly(
  data,
  time_lim = NULL,
  baseline = NULL,
  legend = TRUE,
  allow_facets = FALSE,
  browse_mode = FALSE,
  ...
)

# S3 method for class 'eeg_stats'
plot_butterfly(
  data,
  time_lim = NULL,
  baseline = NULL,
  legend = TRUE,
  allow_facets = FALSE,
  browse_mode = FALSE,
  quantity = "statistic",
  ...
)

# S3 method for class 'eeg_lm'
plot_butterfly(
  data,
  time_lim = NULL,
  baseline = NULL,
  legend = TRUE,
  allow_facets = FALSE,
  browse_mode = FALSE,
  quantity = "coefficients",
  ...
)
```

## Arguments

- data:

  EEG dataset. Should have multiple timepoints.

- ...:

  Other parameters passed to `plot_butterfly`

- time_lim:

  Character vector. Numbers in whatever time unit is used specifying
  beginning and end of time-range to plot. e.g. c(-.1,.3)

- baseline:

  Character vector. Times to use as a baseline. Takes the mean over the
  specified period and subtracts. e.g. c(-.1, 0)

- legend:

  Include plot legend. Defaults to TRUE.

- continuous:

  Is the data continuous or not (I.e. epoched)

- browse_mode:

  Applies custom theme used with
  [`browse_data()`](https://craddm.github.io/eegUtils/reference/browse_data.md).

- allow_facets:

  Allow use of ggplot2 facetting. See note below. Defaults to FALSE.

- quantity:

  Which aspect of the linear model you want to be plotted. only applies
  to `eeg_lm` objects

## Value

ggplot2 object showing ERPs for all electrodes overlaid on a single
plot.

## Methods (by class)

- `plot_butterfly(default)`: Default `plot_butterfly` method for
  data.frames, `eeg_data`

- `plot_butterfly(eeg_evoked)`: Plot butterfly for `eeg_evoked` objects

- `plot_butterfly(eeg_data)`: Create butterfly plot for `eeg_data`
  objects

- `plot_butterfly(eeg_epochs)`: Create butterfly plot for `eeg_epochs`
  objects

- `plot_butterfly(eeg_stats)`: Create butterfly plot for `eeg_stats`
  objects

- `plot_butterfly(eeg_lm)`: Create butterfly plot for `eeg_lm` objects

## Notes on ggplot2 facetting

In order for ggplot2 facetting to work, the data has to be plotted using
[`stat_summary()`](https://ggplot2.tidyverse.org/reference/stat_summary.html)
rather than
[`geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html),
so that the plots can still be made when not all categorical variables
are reflected in the facets. e.g. if there are two variables with two
levels each, but you want to average over one of those variables,
[`stat_summary()`](https://ggplot2.tidyverse.org/reference/stat_summary.html)
is required. However,
[`stat_summary()`](https://ggplot2.tidyverse.org/reference/stat_summary.html)
can be extremely slow with a large number of timepoints, electrodes, and
epochs.

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
plot_butterfly(demo_epochs)
#> Creating epochs based on combinations of variables: epoch_label participant_id 

plot_butterfly(demo_epochs,
time_lim = c(-.1, .4),
legend = FALSE)
#> Creating epochs based on combinations of variables: epoch_label participant_id 
```
