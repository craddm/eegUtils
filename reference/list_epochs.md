# List epochs

List trigger types and any labels found in an `eeg_epochs` object.

## Usage

``` r
list_epochs(data, ...)

# S3 method for class 'eeg_epochs'
list_epochs(data, ...)

# S3 method for class 'eeg_ICA'
list_epochs(data, ...)
```

## Arguments

- data:

  An object of class `eeg_epochs`

- ...:

  Additional arguments

## Methods (by class)

- `list_epochs(eeg_epochs)`: List epochs and associated events from
  `eeg_epochs` objects

- `list_epochs(eeg_ICA)`: List epochs and associated events from
  `eeg_ICA` objects

## See also

[`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)
and
[`list_events()`](https://craddm.github.io/eegUtils/reference/list_events.md)

Other event handlers:
[`events()`](https://craddm.github.io/eegUtils/reference/events.md),
[`list_events()`](https://craddm.github.io/eegUtils/reference/list_events.md),
[`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)

## Author

Matt Craddock <matt@mattcraddock.com>
