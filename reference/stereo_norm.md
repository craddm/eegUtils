# Stereographic electrode projection

Produce a set of x and y coordinates for plotting from 3D Cartesian
coordinates. This is a stereographic projection of the 3D coordinates,
which compensates for the distance of the electrode from the projecting
point and flattens out the scalp.

## Usage

``` r
stereo_norm(chan_info)
```

## Arguments

- chan_info:

  Channel information from an eegUtils objects

## Value

A data.frame with x and y columns indictating electrode positions in
degrees
