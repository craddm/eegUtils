# Create windowed-sinc filter kernel

Supplied with the length of the filter in samples, the cutoff frequency,
and the windowing function, returns a windowed-sinc filter kernel for
use in FIR filtering.

## Usage

``` r
filt_kernel(filt_order, cut_off, w)
```

## Arguments

- filt_order:

  Length of the filter in samples

- cut_off:

  Cutoff frequency (normalized as a fraction of sampling rate)

- w:

  Window
