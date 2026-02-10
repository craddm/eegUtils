# Function to create an S3 object of class `eeg_epochs`.

Function to create an S3 object of class `eeg_epochs`.

## Usage

``` r
eeg_epochs(
  data,
  srate,
  events = NULL,
  chan_info = NULL,
  timings = NULL,
  reference = NULL,
  epochs = NULL
)
```

## Arguments

- data:

  Raw data - signals from electrodes/channels.

- srate:

  Sampling rate in Hz.

- events:

  Event table

- chan_info:

  String of character names for electrodes.

- timings:

  Timing information - samples and sample /sampling rate.

- reference:

  Reference channel information, including names of reference channels,
  excluded channels etc.

- epochs:

  Information about the epochs contained in the data.

## Author

Matt Craddock <matt@mattcraddock.com>
