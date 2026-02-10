# Select electrodes from a given dataset

This is a generic function for selection of electrodes from an EEG
dataset.

## Usage

``` r
select_elecs(data, ...)

# Default S3 method
select_elecs(data, electrode = NULL, keep = TRUE, ...)

# S3 method for class 'eeg_data'
select_elecs(data, electrode, keep = TRUE, df_out = FALSE, ...)

# S3 method for class 'eeg_evoked'
select_elecs(data, electrode = NULL, keep = TRUE, df_out = FALSE, ...)

# S3 method for class 'eeg_ICA'
select_elecs(data, component, keep = TRUE, df_out = FALSE, ...)

# S3 method for class 'eeg_tfr'
select_elecs(data, electrode, keep = TRUE, df_out = FALSE, ...)
```

## Arguments

- data:

  An EEG dataset.

- ...:

  Arguments used with related methods

- electrode:

  A character vector of electrode labels for selection or removal.

- keep:

  Defaults to TRUE. Set to false to *remove* the selected electrodes.

- df_out:

  Defaults to FALSE. Set to TRUE to return a dataframe rather than an
  `eeg_data` object.

- component:

  Component to select

## Value

Data frame with only data from the chosen electrodes

`eeg_data` object with selected electrodes removed/kept.

## Methods (by class)

- `select_elecs(default)`: Select electrodes from a generic data frame.

- `select_elecs(eeg_data)`: Select electrodes from a `eeg_data` object.

- `select_elecs(eeg_evoked)`: Select electrode from an eeg_evoked object

- `select_elecs(eeg_ICA)`: Select components from `eeg_ICA` objects.

- `select_elecs(eeg_tfr)`: Select electrodes from `eeg_tfr` objects.

## See also

[`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)
and
[`select_epochs()`](https://craddm.github.io/eegUtils/reference/select_epochs.md)

Other Data selection functions:
[`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
names(demo_epochs$signals)
#>  [1] "A5"  "A13" "A21" "A29" "A31" "B5"  "B6"  "B8"  "B16" "B18" "B26"
keep_A5 <- select_elecs(demo_epochs, electrode = "A5")
remove_A5 <- select_elecs(demo_epochs, electrode = "A5", keep = FALSE)

select_elecs(demo_epochs, c("A21", "A29"))
#> Epoched EEG data
#> 
#> Number of channels   : 2 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A21 A29 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```
