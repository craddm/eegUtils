# Channel statistics

Channel statistics

## Usage

``` r
channel_stats(data, ...)

# S3 method for class 'eeg_data'
channel_stats(data, ...)
```

## Arguments

- data:

  An `eeg_data` or `eeg_epochs` object.

- ...:

  Other parameters passed to the functions.

## Value

A data frame with statistics for each channel.

## Methods (by class)

- `channel_stats(eeg_data)`: Calculate channel statistics for `eeg_data`
  objects.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
channel_stats(demo_epochs)
#>     electrode      means      sds variance  kurtosis   minmax
#> A5         A5 -2.0201400 5.872054 34.48102 0.8110859 55.99172
#> A13       A13 -0.7884497 3.922905 15.38918 0.8728392 39.18374
#> A21       A21  3.2864957 6.088564 37.07061 1.4236771 68.28431
#> A29       A29  4.1693878 7.781989 60.55936 0.1267196 57.26694
#> A31       A31  1.0943630 5.607857 31.44806 0.8183168 47.41699
#> B5         B5 -2.1856578 7.477209 55.90865 1.0617711 73.42224
#> B6         B6 -2.8632122 5.963742 35.56622 0.3695797 46.77692
#> B8         B8 -1.6502435 6.782543 46.00289 0.9660805 66.80369
#> B16       B16 -1.4560545 5.200351 27.04365 2.0987794 51.42989
#> B18       B18 -0.4696431 5.460927 29.82173 2.3864875 63.04594
#> B26       B26  2.8831543 6.146965 37.78518 0.7989692 56.77661
```
