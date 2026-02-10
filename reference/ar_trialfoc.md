# Detect high trial focality of ICA components

Detect components that load heavily on a small number of trials. Looks
for components that have one particular trial that has a particularly
high z-score.

## Usage

``` r
ar_trialfoc(data, plot = TRUE, threshold = NULL, verbose = TRUE)
```

## Arguments

- data:

  `eeg_ICA` object

- plot:

  Produce plot showing max z-scores and threshold for all ICA
  components.

- threshold:

  Specify a threshold (z-score) for high focality. NULL estimates the
  threshold automatically.

- verbose:

  Print informative messages.

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
ar_trialfoc(demo_sobi)
#> Estimated trial focality threshold (z): 6.41

#> Components with high trial focality: Comp009 
#> [1] "Comp009"
```
