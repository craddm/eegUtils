# Interpolate channels

Interpolate channels

## Usage

``` r
interp_chans(.data, bad_chans, missing_coords = FALSE, weights)
```

## Arguments

- .data:

  Channel data containing all data

- bad_chans:

  Vector of names of bad channels

- missing_coords:

  Logical vector indicating any channels in the data that had no
  associated coordinates

- weights:

  Spherical spline weights for interpolation
