# Object creator for eeg_tfr objects.

Object creator for eeg_tfr objects.

## Usage

``` r
eeg_tfr(
  data,
  srate,
  events,
  chan_info = NULL,
  reference,
  timings = NULL,
  freq_info,
  dimensions,
  epochs = NULL
)
```

## Arguments

- data:

  TFR transformed data

- srate:

  Sampling rate in Hz. A numeric value.

- events:

  Event tables

- chan_info:

  Standard channel information.

- reference:

  Reference information

- timings:

  Timing information.

- freq_info:

  Frequencies and other useful information

- dimensions:

  List of which dimension is which

- epochs:

  Epoch information

## Author

Matt Craddock <matt@mattcraddock.com>
