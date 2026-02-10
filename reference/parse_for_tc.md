# Parse data for timecourses

Internal command for parsing various data structures into a suitable
format for `tc_plot`

## Usage

``` r
parse_for_tc(
  data,
  time_lim,
  electrode,
  baseline,
  add_CI,
  facets,
  mapping,
  colour
)
```

## Arguments

- data:

  data to be parsed

- time_lim:

  time limits to be returned.

- electrode:

  electrodes to be selected

- baseline:

  baseline times to be average and subtracted

- add_CI:

  Logical for whether CIS are required

- facets:

  A RHS-only formula for use with
  [`ggplot2::facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)

- mapping:

  A `ggplot2`
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html) call with
  axis mappings

- colour:

  A character vector indicating which variable to use for colour.
