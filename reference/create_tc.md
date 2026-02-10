# Internal function for creation of timecourse plots

Internal function for creation of timecourse plots

## Usage

``` r
create_tc(
  data,
  add_CI,
  colour,
  quantity = amplitude,
  mapping = NULL,
  facets = NULL
)
```

## Arguments

- data:

  A data frame to be plotted

- add_CI:

  whether to add confidence intervals

- colour:

  whether to use colour

- quantity:

  The name of the column/quantity to plot

- mapping:

  A ggplot2 `aes()` mapping.
