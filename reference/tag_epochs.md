# Tag epochs with labels

Tag epochs with labels indicating details such as experimental
condition, based on the occurrence of event triggers from the events()
structure. This adds a new column to the epochs structure in an
`eeg_epochs` object.

## Usage

``` r
tag_epochs(.data, ...)

# Default S3 method
tag_epochs(.data, ...)

# S3 method for class 'eeg_epochs'
tag_epochs(.data, event_type = NULL, event_label = NULL, ...)
```

## Arguments

- .data:

  An `eegUtils` object

- ...:

  Additional arguments.

- event_type:

  Label epochs according to specific event_types (typically a trigger)

- event_label:

  Label epochs according to specific event_labels

## Methods (by class)

- `tag_epochs(eeg_epochs)`: Tag epochs in an `eeg_epochs` object.
