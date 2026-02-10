# Modify channel information

Get or set the contents of the channel information inside `eegUtils`
objects.

## Usage

``` r
channels(.data)

channels(.data) <- value
```

## Arguments

- .data:

  `eegUtils` object to view

- value:

  Value to replace `chan_info` structure with.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
channels(demo_epochs)
#> # A tibble: 11 × 9
#>    electrode radius theta   phi cart_x cart_y cart_z     x     y
#>    <chr>      <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
#>  1 A5             1   -60   -51  -46.3   57.2  42.5  -37.8  46.6
#>  2 A13            1   -46     0  -61.1    0    59.0  -46     0  
#>  3 A21            1   -60    51  -46.3  -57.2  42.5  -37.8 -46.6
#>  4 A29            1    92   -90    0    -85.0  -2.97   0   -92  
#>  5 A31            1    46   -90    0    -61.1  59.0    0   -46  
#>  6 B5             1    69    90    0     79.4  30.5    0    69  
#>  7 B6             1    46    90    0     61.1  59.0    0    46  
#>  8 B8             1    60    51   46.3   57.2  42.5   37.8  46.6
#>  9 B16            1     0     0    0      0    85      0     0  
#> 10 B18            1    46     0   61.1    0    59.0   46     0  
#> 11 B26            1    60   -51   46.3  -57.2  42.5   37.8 -46.6
```
