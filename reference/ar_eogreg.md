# Remove EOG using regression

Calculates and removes the contribution of eye movements to the EEG
signal using least-squares regression. Specifically, it generate
regression weights based on EOG channels that are used to estimate how
much activity eye movements are responsible for across all channels.

## Usage

``` r
ar_eogreg(data, heog, veog, bipolarize = TRUE)

# S3 method for class 'eeg_data'
ar_eogreg(data, heog, veog, bipolarize = TRUE)

# S3 method for class 'eeg_epochs'
ar_eogreg(data, heog, veog, bipolarize = TRUE)
```

## Arguments

- data:

  Data to regress - `eeg_data` or `eeg_epochs`

- heog:

  Horizontal EOG channel labels

- veog:

  Vertical EOG channel labels

- bipolarize:

  Bipolarize the EOG channels. Only works when four channels are
  supplied (2 HEOG and 2 VEOG).

## Value

An `eeg_data` or `eeg_epochs` object with corrections applied.

## Author

Matt Craddock, <matt@mattcraddock.com>
