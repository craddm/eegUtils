# Function to create an S3 object of class `eeg_ICA`.

Function to create an S3 object of class `eeg_ICA`.

## Usage

``` r
eeg_ICA(
  mixing_matrix,
  unmixing_matrix,
  signals,
  timings,
  events,
  chan_info,
  srate,
  epochs,
  algorithm,
  contents
)
```

## Arguments

- mixing_matrix:

  ICA mixing matrix

- unmixing_matrix:

  ICA unmixing matrix

- signals:

  ICA component timecourses.

- timings:

  Unique timepoints remaining in the data.

- events:

  event table

- chan_info:

  A data frame containing electrode labels and coordinates.

- srate:

  Sampling rate in Hz. A numeric value.

- epochs:

  A data frame containing meta-information about the epochs contained in
  the data, such as participant ID label and condition labels for
  epochs.

- algorithm:

  The method used to calculate the decomposition.

- contents:

  Whether the object contains only the mixing and unmixing matrices
  ("weights") or source timecourses as well ("full")

## Value

An object of class `eeg_ICA`.

## Author

Matt Craddock <matt@mattcraddock.com>
