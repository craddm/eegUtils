# Function to create an object of class eeg_psd

Function to create an object of class eeg_psd

## Usage

``` r
eeg_psd(data, srate, chan_info = NULL, timings = NULL, freqs, epochs)
```

## Arguments

- data:

  PSD transformed data

- srate:

  Sampling rate in Hz.

- chan_info:

  String of character names for electrodes.

- timings:

  Timing information - samples and sample /samplirng rate.

- freqs:

  vector of frequencies

- epochs:

  Epoch information

## Author

Matt Craddock <matt@mattcraddock.com>
