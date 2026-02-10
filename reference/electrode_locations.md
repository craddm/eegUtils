# Get standard electrode locations

Joins standard electrode locations to EEG data from `eegUtils` internal
data.

## Usage

``` r
electrode_locations(data, ...)

# S3 method for class 'data.frame'
electrode_locations(
  data,
  electrode = "electrode",
  drop = FALSE,
  montage = NULL,
  overwrite = TRUE,
  ...
)

# S3 method for class 'eeg_data'
electrode_locations(data, drop = FALSE, montage = NULL, overwrite = FALSE, ...)

# S3 method for class 'eeg_epochs'
electrode_locations(data, drop = FALSE, montage = NULL, overwrite = FALSE, ...)

# S3 method for class 'eeg_tfr'
electrode_locations(data, drop = FALSE, montage = NULL, overwrite = FALSE, ...)
```

## Arguments

- data:

  An EEG dataset.

- ...:

  Passed to S3 methods.

- electrode:

  The column name containing electrode names in data. (Defaults to
  "electrode").

- drop:

  Should electrodes in `data` for which default locations are not
  available be removed? (Defaults to FALSE).

- montage:

  Name of an existing montage set. Defaults to NULL.

- overwrite:

  Overwrite existing channel info. Defaults to FALSE.

## Value

A tibble (or data.frame), or ggplot2 object if `plot = TRUE`.

## Details

The standard locations are from the 10-05 system derived by Oostenveld &
Praamstra (2001). In addition, there are multiple specific montages for
BioSemi systems included. These can be added using the montage
parameter: "biosemi16", "biosemi32", biosemi64", "biosemi64alpha",
"biosemi128", "biosemi160", "biosemi256"

## Methods (by class)

- `electrode_locations(data.frame)`: Adds standard locations to a data
  frame in long format

- `electrode_locations(eeg_data)`: Adds standard locations to the
  chan_info field of an eeg_data object.

- `electrode_locations(eeg_epochs)`: Adds standard locations to the
  chan_info field of an `eeg_data` object.

- `electrode_locations(eeg_tfr)`: Adds standard locations to the
  chan_info field of an `eeg_tfr` object.

## References

Oostenveld, R. & Praamstra, P. (2001). The five percent electrode system
for high-resolution EEG and ERP measurements. Clinical Neurophysiology,
112, 4, 713-719

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
electrode_locations(demo_epochs, overwrite = TRUE, montage = "biosemi64alpha")
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```
