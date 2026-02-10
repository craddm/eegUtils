# Select epochs

This is a generic function for selecting epochs from an epoched data
set.

## Usage

``` r
select_epochs(data, ...)

# Default S3 method
select_epochs(data, ...)

# S3 method for class 'eeg_epochs'
select_epochs(
  data,
  epoch_events = NULL,
  epoch_no = NULL,
  keep = TRUE,
  df_out = FALSE,
  ...
)

# S3 method for class 'eeg_ICA'
select_epochs(
  data,
  epoch_events = NULL,
  epoch_no = NULL,
  keep = TRUE,
  df_out = FALSE,
  ...
)

# S3 method for class 'eeg_tfr'
select_epochs(
  data,
  epoch_events = NULL,
  epoch_no = NULL,
  keep = TRUE,
  df_out = FALSE,
  ...
)
```

## Arguments

- data:

  `eeg_epochs` object from which to select epochs.

- ...:

  Parameters passed to specific methods

- epoch_events:

  Select epochs containing any of the specified events. Can be numeric
  or a character string. Will override any epoch_no input.

- epoch_no:

  Select epochs by epoch number.

- keep:

  Defaults to TRUE, meaning select the specified epochs. Set to FALSE to
  remove specified epochs.

- df_out:

  Output a data.frame instead of an eeg_data object.

## Methods (by class)

- `select_epochs(default)`: Select from generic object

- `select_epochs(eeg_epochs)`: Selection of epochs from `eeg_epochs`
  objects.

- `select_epochs(eeg_ICA)`: Selection of epochs from `eeg_ICA` objects.

- `select_epochs(eeg_tfr)`: Selection of epochs from `eeg_tfr` objects.

## See also

[`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)
and
[`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md)

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
select_epochs(demo_epochs, epoch_no = 1:5)
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 5 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
demo_ica <- run_ICA(demo_epochs, pca = 10)
#> Reducing data to 10 dimensions using PCA.
#> Running SOBI ICA.
select_epochs(demo_ica, epoch_no = 1:5)
#> Epoched ICA decomposition
#> 
#> Number of components : 10 
#> Number of epochs : 5 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Sampling rate        : 128  Hz
```
