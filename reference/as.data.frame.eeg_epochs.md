# Convert `eeg_epochs` object to data.frame

Convert an `eeg_epochs` object to a data.frame for use with whatever
packages you desire.

## Usage

``` r
# S3 method for class 'eeg_epochs'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  long = FALSE,
  events = FALSE,
  cond_labels,
  coords = TRUE,
  ...
)
```

## Arguments

- x:

  Object of class `eeg_epochs`

- row.names:

  Kept for compatability with S3 generic, ignored.

- optional:

  Kept for compatability with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE.

- events:

  Include events in output. Defaults to FALSE. Currently ignored.

- cond_labels:

  Add column tagging epochs with events that have matching labels.
  Deprecated. Metainfo from the epochs structure is now added
  automatically.

- coords:

  Include electrode coordinates in output. Ignored if long == FALSE.

- ...:

  arguments for other as.data.frame commands

## Author

Matt Craddock <matt@mattcraddock.com>
