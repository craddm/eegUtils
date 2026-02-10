# Baseline correction

Used to correct data using the mean of a specified time period. For
time-domain data, this will subtract the mean from all data. For
`eeg_tfr` objects, a variety of methods are available, including
subtraction, and conversion to "dB" change. With a data frame, it will
search for "electrode" and "epoch" columns, and groups on these when
found. An electrode column is always required; an epoch column is not.
Note that baseline correction is always applied on single-trial basis.
For baseline correction based on subtraction, this makes no difference
compared to averaging first and then baseline correcting, but for
divisive measures used with time-frequency data, this distinction can be
very important, and can lead to counterintuitive results.

## Usage

``` r
rm_baseline(data, time_lim = NULL, ...)

# S3 method for class 'eeg_data'
rm_baseline(data, time_lim = NULL, verbose = TRUE, ...)

# S3 method for class 'eeg_epochs'
rm_baseline(data, time_lim = NULL, verbose = TRUE, ...)

# S3 method for class 'data.frame'
rm_baseline(data, time_lim = NULL, verbose = TRUE, ...)

# S3 method for class 'eeg_tfr'
rm_baseline(data, time_lim = NULL, type = "divide", verbose = TRUE, ...)

# S3 method for class 'eeg_evoked'
rm_baseline(data, time_lim = NULL, verbose = TRUE, ...)
```

## Arguments

- data:

  Data to be baseline corrected.

- time_lim:

  Numeric character vector (e.g. time_lim \<- c(-.1, 0)) defining the
  time period to use as a baseline. If the value is NULL, it uses the
  mean of the whole of each epoch if the data is epoched, or the channel
  mean if the data is continuous.

- ...:

  other parameters to be passed to functions

- verbose:

  Defaults to TRUE. Output descriptive messages to console.

- type:

  Type of baseline correction to apply. Options are ("divide", "ratio",
  "absolute", "db", and "pc")

## Value

An `eegUtils` object or a `data.frame`, depending on the input.

## Methods (by class)

- `rm_baseline(eeg_data)`: remove baseline from continuous `eeg_data`

- `rm_baseline(eeg_epochs)`: Remove baseline from `eeg_epochs`

- `rm_baseline(data.frame)`: Legacy method for data.frames

- `rm_baseline(eeg_tfr)`: Method for `eeg_tfr` objects

- `rm_baseline(eeg_evoked)`: Method for `eeg_evoked` objects

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
rm_baseline(demo_epochs)
#> Removing channel means per epoch...
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
rm_baseline(demo_epochs, c(-.1, 0))
#> Baseline: -0.1 - 0s
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```
