# Orthographic electrode projection

Produce a set of x and y coordinates for plotting from 3D Cartesian
coordinates. This is an orthographic projection of the 3D coordinates,
resulting in bunching up of electrodes at the further reaches of the
head.

## Usage

``` r
ortho_norm(chan_info)
```

## Arguments

- chan_info:

  Channel information from an eegUtils objects

## Value

A data.frame with x and y columns indictating electrode positions in mm
