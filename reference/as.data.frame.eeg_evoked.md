# Convert `eeg_evoked` object to data frame

Convert `eeg_evoked` object to data frame

## Usage

``` r
# S3 method for class 'eeg_evoked'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  long = FALSE,
  coords = TRUE,
  ...
)
```

## Arguments

- x:

  Object of class `eeg_evoked`

- row.names:

  Kept for compatability with S3 generic, ignored.

- optional:

  Kept for compatability with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE

- coords:

  include electrode coordinates in output. Ignored if long = FALSE.

- ...:

  arguments for other as.data.frame commands

## Author

Matt Craddock <matt@mattcraddock.com>
