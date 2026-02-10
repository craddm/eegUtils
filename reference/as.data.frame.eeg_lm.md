# Convert `eeg_lm` to `data.frame`

Convert an object of class `eeg_data` into a standard `data.frame`.

## Usage

``` r
# S3 method for class 'eeg_lm'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  long = FALSE,
  coords = TRUE,
  values = c("coefficients", "std_err", "t_stats", "r_sq"),
  ...
)
```

## Arguments

- x:

  Object of class `eeg_lm`

- row.names:

  Kept for compatibility with S3 generic, ignored.

- optional:

  Kept for compatibility with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE.

- coords:

  Include electrode coordinates in output. Only possible when long =
  TRUE.

- values:

  Defaults to "coefficients", returning fitted coefficient values.

- ...:

  arguments for other as.data.frame commands

## Author

Matt Craddock <matt@mattcraddock.com>
