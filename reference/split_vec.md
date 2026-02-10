# Segment data

Split data into segments for Welch PSD. Any leftover data is discarded
(i.e. if seg_length is 256 and signal length is 400, only 1 segment is
returned)

## Usage

``` r
split_vec(vec, seg_length, overlap, detrend = "mean")
```

## Arguments

- vec:

  Data vector to be split up into segments.

- seg_length:

  Length of segments to be FFT'd (in samples).

- overlap:

  Overlap between segments (in samples).

- detrend:

  Detrend segments. Defaults to "mean" - removes mean from each segment.
  Anything else turns off detrending.

## Author

Matt Craddock <matt@mattcraddock.com>
