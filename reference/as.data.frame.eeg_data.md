# Convert `eeg_data` to `data.frame`

Convert an object of class `eeg_data` into a standard data.frame.

## Usage

``` r
# S3 method for class 'eeg_data'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  long = FALSE,
  events = FALSE,
  coords = FALSE,
  ...
)
```

## Arguments

- x:

  Object of class `eeg_data`

- row.names:

  Kept for compatibility with S3 generic, ignored.

- optional:

  Kept for compatibility with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE

- events:

  Include events in output.

- coords:

  Include electrode coordinates in output. Only possible when long =
  TRUE.

- ...:

  arguments for other as.data.frame commands

## Author

Matt Craddock <matt@mattcraddock.com>
