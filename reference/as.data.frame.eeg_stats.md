# Convert `eeg_stats` objects to data frames

Convert `eeg_stats` objects to data frames

## Usage

``` r
# S3 method for class 'eeg_stats'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  long = FALSE,
  coords = FALSE,
  ...
)
```

## Arguments

- x:

  Object of class `eeg_stats`

- row.names:

  Kept for compatability with S3 generic, ignored.

- optional:

  Kept for compatability with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE.

- coords:

  Include electrode coordinates in output (ignored if long = FALSE)

- ...:

  arguments for other as.data.frame commands

## Author

Matt Craddock <matt@mattcraddock.com>
