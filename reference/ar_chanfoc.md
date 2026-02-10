# Detect high channel focality of ICA components

Detect components that load heavily on a single channel. Looks for
components that have one particular channel that has a particularly high
z-score.

## Usage

``` r
ar_chanfoc(
  data,
  plot = TRUE,
  threshold = NULL,
  verbose = TRUE,
  measure = "max",
  ...
)
```

## Arguments

- data:

  An `eeg_ICA` object

- plot:

  Produce plot showing max z-scores and threshold for all ICA
  components.

- threshold:

  Specify a threshold for high focality. NULL estimates the threshold
  automatically.

- verbose:

  Print informative messages.

- measure:

  Use maximum "max" or "kurtosis".

- ...:

  additional parameters

## Value

A character vector of component names that break the threshold.

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
ar_chanfoc(demo_sobi)
#> Estimated threshold: 2.95

#> Components with high channel focality:  
#> character(0)
```
