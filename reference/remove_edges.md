# Remove convolution edges

Create a matrix indicating which timepoints likely suffer from edge
effects. Returns a time by frequency matrix with NA

## Usage

``` r
remove_edges(sigtime, sigma_t)
```

## Arguments

- sigtime:

  timepoints in the signal

- sigma_t:

  standard deviations of the morlet wavelets
