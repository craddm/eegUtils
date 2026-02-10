# Tag events

Give trigger events meaningful labels. Existing labels will be
overwritten. Use hierarchical labelling to tag an event with multiple
labels: separate labels with a "/" symbol. (e.g. "cond1" for a trigger
that belongs to one condition, "cond1/cond2" for a trigger that could
belong to more than one condition).

## Usage

``` r
tag_events(data, ...)

# S3 method for class 'eeg_data'
tag_events(data, trigs, event_label, ...)

# S3 method for class 'eeg_epochs'
tag_events(data, trigs, event_label, ...)
```

## Arguments

- data:

  An object of class `eeg_data` or `eeg_epochs`

- ...:

  Parameters passed to S3 methods

- trigs:

  Character vector of trigger numbers

- event_label:

  Labels for the events.

## Methods (by class)

- `tag_events(eeg_data)`: Tag events in an `eeg_data` object

- `tag_events(eeg_epochs)`: Tag events in an epoched dataset

## See also

Other event handlers:
[`events()`](https://craddm.github.io/eegUtils/reference/events.md),
[`list_epochs()`](https://craddm.github.io/eegUtils/reference/list_epochs.md),
[`list_events()`](https://craddm.github.io/eegUtils/reference/list_events.md)

## Author

Matt Craddock <matt@mattcraddock.com>
