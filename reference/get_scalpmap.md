# Calculate an interpolated scalpmap

Calculate and return an interpolated scalpmap from a `data.frame` or
`eeg_epochs` object. This provides the basis for a raster plot, as used
to create topographical plots.

## Usage

``` r
get_scalpmap(data, ...)

# S3 method for class 'data.frame'
get_scalpmap(
  data,
  method = "biharmonic",
  grid_res = 100,
  interp_limit = "skirt",
  quantity = "amplitude",
  facets = NULL,
  r = 95,
  k = -1,
  ...
)

# S3 method for class 'eeg_epochs'
get_scalpmap(
  data,
  method = "biharmonic",
  grid_res = 100,
  interp_limit = "skirt",
  quantity = "amplitude",
  facets = NULL,
  r = 95,
  ...
)
```

## Arguments

- data:

  An object of class 'data.frame\`

- ...:

  Other arguments

- method:

  "biharmonic" or "gam" smooth

- grid_res:

  Grid resolution

- interp_limit:

  Interpolate up to the "head" or add a "skirt"

- quantity:

  Quantity to be plotted. Defaults to "amplitude".

- facets:

  Any facets you plan to use

- r:

  Size of headshape

- k:

  Degrees of freedom used for spline when using `method = gam`. Defaults
  to -1, which attempts to automatically determine k.

## Value

A `tibble` with the columns `x`, `y`, `fill`

## Methods (by class)

- `get_scalpmap(data.frame)`: Get scalpmap from a `data.frame`

- `get_scalpmap(eeg_epochs)`: Get scalpmap from an `eeg_epochs` object.

## Examples

``` r
head(get_scalpmap(demo_epochs))
#> Joining with `by = join_by(electrode)`
#> # A tibble: 6 × 3
#>       x     y  fill
#>   <dbl> <dbl> <dbl>
#> 1 -113. -9.60  5.77
#> 2 -113. -5.76  5.54
#> 3 -113. -1.92  5.33
#> 4 -113.  1.92  5.14
#> 5 -113.  5.76  4.96
#> 6 -113.  9.60  4.81
```
