# Create epochs from EEG data

Creates epochs around specified event triggers. Requires data of class
`eeg_data`. Where multiple events are specified, epochs will be created
around each event.

## Usage

``` r
epoch_data(data, ...)

# S3 method for class 'eeg_data'
epoch_data(
  data,
  events,
  time_lim = c(-1, 1),
  baseline = "none",
  epoch_labels = NULL,
  ...
)
```

## Arguments

- data:

  Continuous data to be epoched.

- ...:

  Parameters passed to functions

- events:

  Character vector of events to epoch around.

- time_lim:

  Time in seconds to form epoch around the events. Defaults to one
  second either side.

- baseline:

  Baseline times to subtract. Can be set to a numeric vector of length
  two to specify a time window to use as a baseline in each epoch (e.g.
  c(-.1, 0)), "none", which will perform no baseline correction, or
  `NULL` to use the mean of the whole epoch. As of v0.6 of `eegUtils`,
  the default is "none".

- epoch_labels:

  Character vector of same length as events which'll be used to label
  the epochs.

## Value

Returns an epoched object of class `eeg_epochs`

## Methods (by class)

- `epoch_data(eeg_data)`: Epoch `eeg_data` objects

## Author

Matt Craddock <matt@mattcraddock.com>
