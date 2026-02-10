# Convert `eeg_tfr` objects to a `data.frame`

Convert `eeg_tfr` objects to a `data.frame`

## Usage

``` r
# S3 method for class 'eeg_tfr'
as.data.frame(x, row.names = NULL, optional = FALSE, long = FALSE, ...)
```

## Arguments

- x:

  Object of class `eeg_tfr`

- row.names:

  Kept for compatability with S3 generic, ignored.

- optional:

  Kept for compatability with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE.

- ...:

  arguments for other `as.data.frame` commands

## Value

A `data.frame` or `tibble`

## Author

Matt Craddock <matt@mattcraddock.com>
