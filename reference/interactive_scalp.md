# Interactive scalp maps

Launches a Shiny Gadget for an interactive version of `erp_scalp`,
allowing clicking of individual electrodes to plot them in a separate
panel. In that panel, they can be averaged over or plotted as individual
electrodes.

## Usage

``` r
interactive_scalp(data, ...)

# S3 method for class 'eeg_epochs'
interactive_scalp(data, colour = NULL, baseline = NULL, montage = NULL, ...)
```

## Arguments

- data:

  An EEG dataset.

- ...:

  Additional arguments

- colour:

  Variable to colour lines by. If no variable is passed, only one line
  is drawn for each electrode.

- baseline:

  Character vector of times to subtract for baseline correction.

- montage:

  Name of an existing montage set. Defaults to NULL; (currently only
  'biosemi64alpha' available other than default 10/20 system)

## Methods (by class)

- `interactive_scalp(eeg_epochs)`: Method for `eeg_epochs` objects

## See also

[`erp_scalp()`](https://craddm.github.io/eegUtils/reference/erp_scalp.md)
for non-interactive plots of ERPs in a scalp-based layout.

Other scalp-based maps:
[`erp_scalp()`](https://craddm.github.io/eegUtils/reference/erp_scalp.md),
[`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)

## Author

Matt Craddock, <matt@mattcraddock.com>
