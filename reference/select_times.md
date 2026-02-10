# Select timerange

Generic function for selecting specific time ranges from a given
dataset. Input can be a dataframe, or an object of class `eeg_data`,
`eeg_epochs`, or `eeg_evoked`. Note this finds the closest times to
those specified, so the time range returned may be slightly longer or
shorter than that requested.

## Usage

``` r
select_times(data, ...)

# Default S3 method
select_times(data, time_lim = NULL, ...)

# S3 method for class 'eeg_data'
select_times(data, time_lim = NULL, df_out = FALSE, ...)

# S3 method for class 'eeg_epochs'
select_times(data, time_lim, df_out = FALSE, ...)

# S3 method for class 'eeg_evoked'
select_times(data, time_lim, df_out = FALSE, ...)

# S3 method for class 'eeg_tfr'
select_times(data, time_lim = NULL, df_out = FALSE, ...)
```

## Arguments

- data:

  Data from which to select

- ...:

  Further arguments passed to or from other methods.

- time_lim:

  A character vector of two numbers indicating the time range to be
  selected e.g. c(min, max), or a list to select specific times.

- df_out:

  Returns a data frame rather than an object of the same type that was
  passed in.

## Value

Data frame with only data from within the specified range.

`eeg_data` object

## Methods (by class)

- `select_times(default)`: Default select times function

- `select_times(eeg_data)`: Select times from an `eeg_data` object

- `select_times(eeg_epochs)`: Select times in `eeg_epochs` objects

- `select_times(eeg_evoked)`: Select times in `eeg_evoked` objects

- `select_times(eeg_tfr)`: Select times from an `eeg_tfr` object

## See also

[`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md)
and
[`select_epochs()`](https://craddm.github.io/eegUtils/reference/select_epochs.md)

Other Data selection functions:
[`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md)

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
## Select timepoints from -.1 to .3
demo_epochs
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
short_epochs <- select_times(demo_epochs, time_lim = c(-.1, .3))
short_epochs
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.096 - 0.295 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
short_epochs <- select_times(demo_epochs, time_lim = list(-.1, .15, .3))
#> Returning closest time points to those requested: -0.104 0.146 0.295 
head(short_epochs$timings)
#> # A tibble: 6 × 3
#>   epoch sample   time
#>   <dbl>  <dbl>  <dbl>
#> 1     1   4076 -0.104
#> 2     1   4204  0.146
#> 3     1   4280  0.295
#> 4     2   6985 -0.104
#> 5     2   7113  0.146
#> 6     2   7189  0.295
```
