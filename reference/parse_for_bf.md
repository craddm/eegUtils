# Parse data for butterfly plots

Internal command for parsing various data structures into a suitable
format for `plot_butterfly`

## Usage

``` r
parse_for_bf(data, time_lim = NULL, baseline = NULL, quantity = "coefficients")
```

## Arguments

- data:

  data to be parsed

- time_lim:

  time limits to be returned.

- baseline:

  baseline times to be average and subtracted
