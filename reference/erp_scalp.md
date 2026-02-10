# Plot event-related potentials using a scalp based layout

Creates a `ggplot2` figure showing the ERP at each electrode, with each
electrode placed according to its location on the scalp.

## Usage

``` r
erp_scalp(data, ...)

# Default S3 method
erp_scalp(
  data,
  electrode = "electrode",
  amplitude = "amplitude",
  time = "time",
  color = NULL,
  colour = NULL,
  size = 0.8,
  linewidth = 0.8,
  baseline = NULL,
  show_guide = TRUE,
  chan_info = NULL,
  montage = NULL,
  ...
)
```

## Arguments

- data:

  An `eeg_epochs` or `eeg_evoked` object

- ...:

  Various arguments passed to specific functions

- electrode:

  Column name containing electrode names in data. Defaults to
  "electrode".

- amplitude:

  Column name containing amplitudes in data. Defaults to "amplitude".

- time:

  Column name containing time in data. Defaults to "time".

- color:

  Alias for `colour`.

- colour:

  Variable to colour lines by. If no variable is passed, only one line
  is drawn for each electrode.

- size:

  Size of the line(s) for the ERPs. Deprecated

- linewidth:

  Size of the line(s) for the ERPs.

- baseline:

  Character vector of times to subtract for baseline correct.

- show_guide:

  Should be a guide showing the scale of the ERP plots be shown.
  Defaults to TRUE.

- chan_info:

  Pass channel information in the standard `eegUtils` format directly.

- montage:

  Name of an existing montage set. Defaults to NULL, which will attempt
  to us locations from the data. If none are found, it will attempt to
  use standard 10-05 locations.

## Value

Returns a `ggplot2` plot object.

## Methods (by class)

- `erp_scalp(default)`: Plot an ERP scalp map

## See also

[`interactive_scalp()`](https://craddm.github.io/eegUtils/reference/interactive_scalp.md)
for interactive plots of ERPs in a scalp-based layout.

Other scalp-based maps:
[`interactive_scalp()`](https://craddm.github.io/eegUtils/reference/interactive_scalp.md),
[`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)

## Author

Matti Vuorre, <mv2521@columbia.edu>

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
erp_scalp(demo_epochs, montage = "biosemi64alpha")

```
