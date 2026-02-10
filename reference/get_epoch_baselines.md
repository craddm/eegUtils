# Get epoch baselines

Gets the baseline values for every epoch separately

## Usage

``` r
get_epoch_baselines(data, time_lim)
```

## Arguments

- data:

  data for which to calculate the baselines

- time_lim:

  time limits of the baseline period. numeric vector of length two,
  c(start, end)

## Value

A numeric matrix of n_epochs x n_channels.
