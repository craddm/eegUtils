# Compute power spectral density

`compute_psd` returns the PSD calculated using Welch's method for every
channel in the data. The output is in microvolts-squared divided by
Hertz - \\\muV^2 / Hz\\. If the object has multiple epochs, it will
perform Welch's FFT separately for each epoch and then average them
afterwards.

## Usage

``` r
compute_psd(data, ...)

# S3 method for class 'eeg_data'
compute_psd(
  data,
  seg_length = NULL,
  noverlap = NULL,
  n_fft = NULL,
  method = "Welch",
  demean = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'eeg_epochs'
compute_psd(
  data,
  seg_length = NULL,
  noverlap = NULL,
  n_fft = 256,
  method = "Welch",
  keep_trials = TRUE,
  demean = TRUE,
  verbose = TRUE,
  ...
)

# S3 method for class 'eeg_evoked'
compute_psd(
  data,
  seg_length = NULL,
  noverlap = NULL,
  n_fft = 256,
  method = "Welch",
  demean = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Data to be plotted. Accepts objects of class `eeg_data`

- ...:

  any further parameters passed to specific methods

- seg_length:

  Length of rolling data segments. Defaults to `n_fft`. Must be \<=
  `n_fft`.

- noverlap:

  Number of (sampling) points of overlap between segments. Must be \<=
  `seg_length`.

- n_fft:

  Length of FFT to be calculated in sampling points. See details.

- method:

  Defaults to "Welch". No other method currently implemented.

- demean:

  Remove channel/epoch means. TRUE by default.

- verbose:

  Print informative messages. TRUE by default.

- keep_trials:

  Include FFT for every trial in output, or average over them if FALSE.

## Value

Currently, a data frame with the PSD for each channel separately.

## Details

Welch's FFT splits the data into multiple segments, calculates the FFT
separately for each segment, and then averages over segments. Each
segment is windowed with a Hanning window to counter spectral leakage.
For epoched data, Welch's FFT is calculated separately for each trial.

The number of sampling points used for the FFT can be specified using
n_fft. n_fft defaults to 256 sampling points for `eeg_epochs` data, or
the minimum of 2048 or the length of the signal for continuous
`eeg_data`.

`seg_length` defaults to be `n_fft`, and must be less than or equal to
it.

`noverlap` specifies the amount of overlap between windows in sampling
points. If NULL, it defaults to 50\\

## Methods (by class)

- `compute_psd(eeg_data)`: Compute PSD for an `eeg_data` object

- `compute_psd(eeg_epochs)`: Compute PSD for an `eeg_epochs` object

- `compute_psd(eeg_evoked)`: Compute PSD for an `eeg_evoked` object

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
out <- compute_psd(demo_epochs)
#> Removing channel means per epoch...
#> Computing Power Spectral Density using Welch's method.
#> FFT length: 256
#> Segment length: 84
#> Overlapping points: 42 (50% overlap)

out <- compute_psd(demo_epochs, n_fft = 256, seg_length = 128)
#> Removing channel means per epoch...
#> Computing Power Spectral Density using Welch's method.
#> FFT length: 256
#> Segment length: 84
#> Overlapping points: 42 (50% overlap)
```
