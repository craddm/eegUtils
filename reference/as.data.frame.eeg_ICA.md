# Convert `eeg_ICA` object to data frame

Convert `eeg_ICA` object to data frame

## Usage

``` r
# S3 method for class 'eeg_ICA'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  long = FALSE,
  cond_labels,
  mixing = FALSE,
  coords = TRUE,
  ...
)
```

## Arguments

- x:

  Object of class `eeg_ICA`

- row.names:

  Kept for compatibility with S3 generic, ignored.

- optional:

  Kept for compatibility with S3 generic, ignored.

- long:

  Convert to long format. Defaults to FALSE

- cond_labels:

  add condition labels to data frame. Deprecated.

- mixing:

  If TRUE, outputs the mixing matrix. If FALSE, outputs source
  activations.

- coords:

  Adds electrode coordinates if TRUE; only if long data and the mixing
  matrix are requested.

- ...:

  arguments for other as.data.frame commands

## Author

Matt Craddock <matt@mattcraddock.com>
