# Referencing

Used to reference the EEG data to a specified electrode or electrodes.
Defaults to average reference. When specific electrodes are used, they
are removed from the data. Meta-data about the referencing scheme is
held in the `eeg_data` structure.

## Usage

``` r
eeg_reference(data, ...)

# Default S3 method
eeg_reference(data, ...)

# S3 method for class 'eeg_data'
eeg_reference(
  data,
  ref_chans = "average",
  exclude = NULL,
  robust = FALSE,
  implicit_ref = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'eeg_epochs'
eeg_reference(
  data,
  ref_chans = "average",
  exclude = NULL,
  robust = FALSE,
  implicit_ref = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Data to re-reference. Primarily meant for use with data of class
  `eeg_data`.

- ...:

  Further parameters to be passed to `eeg_reference`

- ref_chans:

  Channels to reference data to. Defaults to "average" i.e. average of
  all electrodes in data. Character vector of channel names or numbers.

- exclude:

  Electrodes to exclude from average reference calculation.

- robust:

  Use median instead of mean; only used for average reference. Defaults
  to FALSE.

- implicit_ref:

  Implicit reference channel - use this to add a channel back that was
  previously used as a reference. E.g. if the LM (left mastoid) channel
  was used in recording and is absent from the data, passing "LM" adds
  an "LM" channel back to the data, populated with zeroes.

- verbose:

  Print informative messages in console. Defaults to TRUE.

## Value

object of class `eeg_data`, re-referenced as requested.

## Methods (by class)

- `eeg_reference(default)`: Default method

- `eeg_reference(eeg_data)`: Rereference objects of class `eeg_data`

- `eeg_reference(eeg_epochs)`: Rereference objects of class `eeg_epochs`

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
# demo_epochs is average referenced by default
demo_epochs
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
# Rereference it but exclude B5 from calculation of the average
eeg_reference(demo_epochs, exclude = "B5")
#> You have used the existing reference channel(s), average again.
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
# Reference data using the median of the reference channels rather than the mean
eeg_reference(demo_epochs, robust = TRUE)
#> You have used the existing reference channel(s), average again.
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 

eeg_reference(demo_spatial)
#> You have used the existing reference channel(s), average again.
#> Epoched EEG data
#> 
#> Number of channels   : 23 
#> Number of epochs : 80 
#> Epoch limits     : -0.301 - 0.488 seconds
#> Electrode names      : Fp1 Fp2 Fz Cz Pz P3 P4 C3 C4 P7 P8 F3 F4 T7 T8 F7 F8 Oz Fpz EXG1 EXG2 EXG3 EXG4 
#> Sampling rate        : 128  Hz
#> Reference        : average 
eeg_reference(demo_spatial, ref_chans = "Fz")
#> Epoched EEG data
#> 
#> Number of channels   : 22 
#> Number of epochs : 80 
#> Epoch limits     : -0.301 - 0.488 seconds
#> Electrode names      : Fp1 Fp2 Cz Pz P3 P4 C3 C4 P7 P8 F3 F4 T7 T8 F7 F8 Oz Fpz EXG1 EXG2 EXG3 EXG4 
#> Sampling rate        : 128  Hz
#> Reference        : Fz 
eeg_reference(demo_spatial, implicit_ref = "LM")
#> You have used the existing reference channel(s), average again.
#> Epoched EEG data
#> 
#> Number of channels   : 24 
#> Number of epochs : 80 
#> Epoch limits     : -0.301 - 0.488 seconds
#> Electrode names      : Fp1 Fp2 Fz Cz Pz P3 P4 C3 C4 P7 P8 F3 F4 T7 T8 F7 F8 Oz Fpz EXG1 EXG2 EXG3 EXG4 LM 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```
