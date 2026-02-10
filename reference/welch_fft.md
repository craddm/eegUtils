# Welch FFT

Internal function for calculating the PSD using Welch's method.
Hardcoded to use Hamming window.

## Usage

``` r
welch_fft(data, seg_length, n_fft, noverlap, n_sig, srate)
```

## Arguments

- data:

  Object to perform FFT on.

- seg_length:

  length of each segment of data.

- n_fft:

  length of FFT.

- noverlap:

  overlap between segments.

- n_sig:

  number of samples total.

- srate:

  Sampling rate of the data.
