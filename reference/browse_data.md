# Browse EEG data

A Shiny gadget for browsing EEG data and ICA decompositions
interactively. With EEG data (epoched or continuous), data can be viewed
as a butterfly plot (all electrodes overlaid) or as individual traces
(electrodes "stacked"). For `eeg_ICA` objects, you will instead be shown
a composite of multiple properties of the decomposition - a topography,
an ERP image, an ERP, and a power spectral density plot from 4-50 Hz.

## Usage

``` r
browse_data(data, ...)

# S3 method for class 'eeg_ICA'
browse_data(data, ...)

# S3 method for class 'eeg_data'
browse_data(data, sig_length = 5, downsample = TRUE, ...)

# S3 method for class 'eeg_epochs'
browse_data(data, sig_length = 5, downsample = FALSE, ...)
```

## Arguments

- data:

  `eeg_data`, `eeg_epochs`, or `eeg_ICA` object to be plotted.

- ...:

  Other parameters passed to browsing functions.

- sig_length:

  Length of signal to be plotted initially (seconds if continuous,
  epochs if epoched).

- downsample:

  Only works on `eeg_data` or `eeg_epochs` objects. Reduces size of data
  by only plotting every 4th point, speeding up plotting considerably.
  Defaults to TRUE for `eeg_data`, FALSE for `eeg_epochs`

## Value

A character vector of component names selected for rejection.

## Methods (by class)

- `browse_data(eeg_ICA)`: View `eeg_ICA` component properties

- `browse_data(eeg_data)`: Browse continuous EEG data.

- `browse_data(eeg_epochs)`: Browse epoched EEG data.

## Author

Matt Craddock <matt@mattcraddock.com>
