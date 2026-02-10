# Detect low autocorrelation of ICA components

Low autocorrelation can be a sign of a poor quality channel or
component. Often these are noisy, poor contact, or heavily contaminated
with muscle noise. Low autocorrelation at a lag of 20ms is often
associated with muscle noise.

## Usage

``` r
ar_acf(data, ...)

# S3 method for class 'eeg_ICA'
ar_acf(data, ms = 20, plot = TRUE, verbose = TRUE, threshold = NULL, ...)
```

## Arguments

- data:

  `eeg_ICA` object

- ...:

  additional parameters

- ms:

  Time lag to check ACF, in milliseconds. Defaults to 20 ms.

- plot:

  Produce plot showing ACF and threshold for all EEG components.

- verbose:

  Print informative messages. Defaults to TRUE.

- threshold:

  Specify a threshold for low ACF. NULL estimates the threshold
  automatically.

## Value

A character vector of component names that break the threshold.

## Methods (by class)

- `ar_acf(eeg_ICA)`: Autocorrelation checker for `eeg_ICA` objects

## References

Chaumon, M., Bishop, D.V., Busch, N.A. (2015). A practical guide to the
selection of independent components of the electroencephalogram for
artifact correction. J Neurosci Methods. Jul 30;250:47-63. doi:
10.1016/j.jneumeth.2015.02.025

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
demo_sobi <- run_ICA(demo_epochs, pca = 10)
#> Reducing data to 10 dimensions using PCA.
#> Running SOBI ICA.
ar_acf(demo_sobi)
#> Estimating autocorrelation at 20ms lag.
#> Estimated ACF threshold: -0.16

#> Subthreshold components:  
#> character(0)
```
