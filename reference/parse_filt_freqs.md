# Parse filter frequency input

Parses the frequencies input by the user, converting them to a fraction
of the sampling rate and setting the filter type (low-pass, high-pass,
band-pass, band-stop) appropriately.

## Usage

``` r
parse_filt_freqs(low_freq, high_freq, srate, method)
```

## Arguments

- low_freq:

  low frequency cutoff (Hz)

- high_freq:

  High frequency cutoff (Hz)

- srate:

  Sampling rate (Hz)

- method:

  "iir" or "fir" method.
