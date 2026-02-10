# Select frequencies

Select specific frequencies from `eeg_tfr` objects. Can be used to
selecting either single frequencies or anything within a range.

## Usage

``` r
select_freqs(data, freq_range)

# S3 method for class 'eeg_tfr'
select_freqs(data, freq_range)
```

## Arguments

- data:

  An `eeg_tfr` object.

- freq_range:

  The range of frequencies to retain. Can be a scale or the lower and
  upper bounds. (e.g. c(5, 30))

## Methods (by class)

- `select_freqs(eeg_tfr)`: Function for selecting specific frequencies
  from `eeg_tfr` objects.

## Examples

``` r
demo_tfr <- compute_tfr(demo_epochs, foi = c(4, 30), n_freq = 10, n_cycles = 5)
#> Computing TFR using Morlet wavelet convolution
#> Output frequencies using linear spacing: 4 6.89 9.78 12.67 15.56 18.44 21.33 24.22 27.11 30
#> Removing channel means per epoch...
#> Returning signal averaged over all trials.
select_freqs(demo_tfr, c(8, 12))
#> Epoched EEG TFR data
#> 
#> Frequency range      :    9.78 
#> Number of channels   :    11 
#> Electrode names      :    A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Number of epochs :    1 
#> Epoch limits     :    -0.197 - 0.451 seconds
#> Sampling rate        :    128  Hz
```
