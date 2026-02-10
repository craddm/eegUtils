# Calculate cycles

A helper function for calculating the appropriate min-max cycles for a
fixed time window/frequency resolution for use with `compute_tfr`. For
some analyses you may wish to keep a fixed frequency resolution across
the range being analysed, which requires using a fixed time window.
`compute_tfr` expects the minimum and maximum number of cycles to be
supplied. Use this function to calculate the equivalent number of cycles
at each frequency.

## Usage

``` r
cycle_calc(time_win, frex)
```

## Arguments

- time_win:

  Time window in seconds.

- frex:

  Frequencies of interest.

## Value

the number of cycles for each frequency of interest

## See also

[`compute_tfr`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)

## Examples

``` r
cycle_calc(.5, seq(3, 30, length.out = 10))
#>  [1]  1.5  3.0  4.5  6.0  7.5  9.0 10.5 12.0 13.5 15.0
no_scale_tfr <- compute_tfr(demo_epochs, foi = c(3, 30),
 n_cycles = range(cycle_calc(0.5, seq(3, 30, length.out = 10))),
  n_freq = 10)
#> Computing TFR using Morlet wavelet convolution
#> Output frequencies using linear spacing: 3 6 9 12 15 18 21 24 27 30
#> Removing channel means per epoch...
#> Returning signal averaged over all trials.
```
