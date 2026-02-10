# Generate filter coefficients

Generate filter coefficients for an IIR or FIR filter.

## Usage

``` r
filter_coefs(method, filt_pars, filter_order, window)
```

## Arguments

- method:

  IIR or FIR

- filt_pars:

  output of parse_filt_freqs

- filter_order:

  order of the filter in samples

- window:

  Ignored for IIR
