# Detect high component correlation with eye channels

Checks the correlation between individual components of an `eeg_ICA`
decomposition and the electrooculogram channels of an `eeg_epochs`
dataset.

## Usage

``` r
ar_eogcor(decomp, data, ...)

# S3 method for class 'eeg_ICA'
ar_eogcor(
  decomp,
  data,
  HEOG,
  VEOG,
  threshold = NULL,
  plot = TRUE,
  bipolarize = TRUE,
  method = c("pearson", "kendall", "spearman"),
  verbose = TRUE,
  ...
)
```

## Arguments

- decomp:

  An `eeg_ica` object

- data:

  The original `eeg_epochs` object from which this decomposition was
  derived

- ...:

  Other parameters

- HEOG:

  Horizontal eye channels

- VEOG:

  Vertical eye channels

- threshold:

  Threshold for correlation (r). Defaults to NULL, automatically
  determining a threshold.

- plot:

  Plot correlation coefficient for all components

- bipolarize:

  Bipolarize the HEOG and VEOG channels?

- method:

  Correlation method. Defaults to Pearson.

- verbose:

  Print informative messages. Defaults to TRUE.

## Value

A character vector of component names that break the threshold.

## Methods (by class)

- `ar_eogcor(eeg_ICA)`: Method for `eeg_ICA` objects.

## References

Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide to the
selection of independent components of the electroencephalogram for
artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
10.1016/j.jneumeth.2015.02.025

## Author

Matt Craddock, <matt@mattcraddock.com>
