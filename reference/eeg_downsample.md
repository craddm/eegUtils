# Downsampling EEG data

Performs low-pass anti-aliasing filtering and downsamples EEG data by a
specified factor. This is a wrapper for `decimate` from the `signal`
package. Note that this will also adjust the event table, moving events
to the nearest time remaining after downsampling

## Usage

``` r
eeg_downsample(data, ...)

# S3 method for class 'eeg_data'
eeg_downsample(data, q, ...)

# S3 method for class 'eeg_epochs'
eeg_downsample(data, q, ...)
```

## Arguments

- data:

  An `eeg_data` object to be downsampled

- ...:

  Parameters passed to functions

- q:

  Integer factor to downsample by

## Methods (by class)

- `eeg_downsample(eeg_data)`: Downsample eeg_data objects

- `eeg_downsample(eeg_epochs)`: Downsample eeg_epochs objects

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
eeg_downsample(demo_epochs, 2)
#> Downsampling from 128Hz to 64Hz.
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.443 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 64  Hz
#> Reference        : average 
```
