# Plot electrode locations

Produces either a 2D plot of the electrode locations or an interactive
plot of electrode locations in 3D space.

## Usage

``` r
plot_electrodes(data, interact = FALSE)

# Default S3 method
plot_electrodes(data, interact = FALSE)

# S3 method for class 'eeg_data'
plot_electrodes(data, interact = FALSE)

# S3 method for class 'eeg_evoked'
plot_electrodes(data, interact = FALSE)

# S3 method for class 'eeg_tfr'
plot_electrodes(data, interact = FALSE)
```

## Arguments

- data:

  Data with associated electrode locations to be plotted.

- interact:

  Choose 2D cartesian layout, or, if set to TRUE, an interactive 3D plot
  of electrode locations. Defaults to FALSE.

## Value

A `ggplot` or `plotly` figure showing the locations of the electrodes

## Methods (by class)

- `plot_electrodes(default)`: generic plot electrodes function

- `plot_electrodes(eeg_data)`: Plot electrodes associated with an
  `eeg_data` object.

- `plot_electrodes(eeg_evoked)`: Plot electrodes associated with an
  `eeg_evoked` object.

- `plot_electrodes(eeg_tfr)`: Plot electrodes associated with an
  `eeg_data` object.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
plot_electrodes(demo_epochs)

```
