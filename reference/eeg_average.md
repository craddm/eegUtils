# Calculate averages (e.g. event-related potentials) for single datasets

This function is used to create an `eeg_evoked` object from
`eeg_epochs`. By default, it will try to keep different conditions in
the data separate using the `epochs` metadata from the object, thus
yielding one average per condition. Alternatively, the user can specify
which averages they want using the `cols` argument.

## Usage

``` r
eeg_average(data, ...)

# Default S3 method
eeg_average(data, ...)

# S3 method for class 'eeg_epochs'
eeg_average(data, cols = NULL, verbose = TRUE, ...)

# S3 method for class 'eeg_evoked'
eeg_average(data, cols = NULL, weighted = TRUE, verbose = TRUE, ...)

# S3 method for class 'eeg_tfr'
eeg_average(data, cols = NULL, weighted = TRUE, verbose = TRUE, ...)
```

## Arguments

- data:

  An `eeg_epochs` of `eeg_tfr` object.

- ...:

  Other arguments passed to the averaging functions

- cols:

  Columns from the `epochs` structure that the average should group on.
  NULL, the default, uses all columns other than the `epoch` column.

- verbose:

  Print informative messages during averaging. Defaults to TRUE

- weighted:

  Produce a weighted average over epochs, which accounts for upstream
  differences in the number of epochs that contribute to each average.

## Value

An object of class `eeg_evoked` if applied to `eeg_epochs`; `eeg_tfr` if
applied to `eeg_tfr`.

## Methods (by class)

- `eeg_average(default)`: Default method for averaging EEG objects

- `eeg_average(eeg_epochs)`: Create evoked data from `eeg_epochs`objects

- `eeg_average(eeg_evoked)`: average an `eeg_evoked` object over epochs.

- `eeg_average(eeg_tfr)`: average an `eeg_tfr` object over epochs.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
eeg_average(demo_spatial)
#> Creating epochs based on combinations of variables: participant_id epoch_labels 
#> Evoked EEG data
#> 
#> Number of channels   :    23 
#> Epoch limits     : -0.301 - 0.488 seconds
#> Electrode names      :    Fp1 Fp2 Fz Cz Pz P3 P4 C3 C4 P7 P8 F3 F4 T7 T8 F7 F8 Oz Fpz EXG1 EXG2 EXG3 EXG4 
#> Sampling rate        :    128  Hz
eeg_average(demo_spatial, cols = "everything")
#> Creating epochs based on combinations of variables: participant_id 
#> Evoked EEG data
#> 
#> Number of channels   :    23 
#> Epoch limits     : -0.301 - 0.488 seconds
#> Electrode names      :    Fp1 Fp2 Fz Cz Pz P3 P4 C3 C4 P7 P8 F3 F4 T7 T8 F7 F8 Oz Fpz EXG1 EXG2 EXG3 EXG4 
#> Sampling rate        :    128  Hz
```
