# Perform Hanning time-frequency analysis

Internal function for performing Hanning wavelet transforms using
convolution in frequency domain

## Usage

``` r
tf_hanning(
  data,
  foi,
  n_freq,
  spacing,
  n_cycles,
  keep_trials,
  output,
  downsample,
  trim_edges = trim_edges,
  demean = TRUE,
  verbose
)
```

## Arguments

- data:

  Data in `eeg_epochs` format.

- foi:

  Frequencies of interest. Scalar or character vector of the lowest and
  highest frequency to resolve.

- n_freq:

  Number of frequencies to be resolved.

- spacing:

  Use linear or log spacing for frequencies.

- n_cycles:

  Number of cycles at each frequency.

- keep_trials:

  Keep single trials or average over them before returning.

- output:

  Sets whether output is power, phase, or fourier coefficients.

- downsample:

  Downsampling factor (integer).

- demean:

  Remove mean before transforming.

- verbose:

  Print informative messages in console.
