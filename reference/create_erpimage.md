# Function for creating an ERP image

Function for creating an ERP image

## Usage

``` r
create_erpimage(data, electrode, smoothing, clim, interpolate = FALSE, na.rm)
```

## Arguments

- data:

  Data frame to be plotted. Requires an amplitude column.

- electrode:

  electrode at which to generate an ERP image.

- smoothing:

  Number of trials to smooth over when generating image

- clim:

  Character vector of min and max values of plotting colour range. e.g.
  c(-5,5). Defaults to min and max.

- interpolate:

  Turn on `geom_raster()` interpolation for smoother images.

- na.rm:

  Remove trials with NA amplitudes after smoothing. Defaults to TRUE.
