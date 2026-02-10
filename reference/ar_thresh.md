# Simple absolute value thresholding

Reject data based on a simple absolute amplitude threshold. This marks
any timepoint from any electrode.

## Usage

``` r
ar_thresh(data, threshold, reject = FALSE)

# S3 method for class 'eeg_data'
ar_thresh(data, threshold, reject = FALSE)

# S3 method for class 'eeg_epochs'
ar_thresh(data, threshold, reject = FALSE)
```

## Arguments

- data:

  An object of class `eeg_data` or `eeg_epochs`.

- threshold:

  In microvolts. If one value is supplied, it will be treated as a +-
  value.

- reject:

  If TRUE, remove marked data immediately, otherwise mark for
  inspection/rejection. Defaults to FALSE.

## Value

An object of class `eeg_data` or `eeg_epochs`

## Methods (by class)

- `ar_thresh(eeg_data)`: Reject data using a simple threshold.

- `ar_thresh(eeg_epochs)`: Reject data using a simple threshold.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
ar_thresh(demo_epochs, c(100))
#> 0 (0%) samples above 100 uV threshold.
#> 0 (0%) samples below -100 uV threshold.
#> No epochs contain samples above threshold.
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```
